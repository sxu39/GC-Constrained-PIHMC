#include "box.h"
#include <algorithm>
using namespace cpihmc;
using namespace std;

cpihmc::box::box(const input &Input):Cell(Input), LatVec(Cell.LatVec * Cell.LatConst), ElecNum(Cell.ElecNum), MassScal(Input.get_mass_scal()), NType(Cell.NType), ElementType(Input.get_element_type()), NAtoms(Cell.NAtoms), Atoms(Cell.Atoms), Elements(Cell.Elements), Masks(Cell.Masks), AtomIndexes(NAtoms), NBead(Input.get_n_bead()), BeadIndexes(Input.get_bead_index()), Temp(Input.get_temp()), KinEng(0.0), PotEng(0.0), TotEng(0.0), QtmKinEng(0.0)
{
    for (index_t i = 0 ; i < NAtoms ; ++i) AtomIndexes[i] = i;
    for (index_t MaskIndex = Masks.size() - 1 ; MaskIndex >= 0 ; --MaskIndex) AtomIndexes.erase(AtomIndexes.begin()+Masks[MaskIndex]);
    for (const auto &Index : BeadIndexes)
    {
        bool_t IsValid = true;
        if (Index >= NAtoms) IsValid = false;
        else for (const auto &MaskIndex : Masks) if (Index == MaskIndex){IsValid = false; break;}
        if (!IsValid) throw invalid_argument("The atom with the index "+to_string(Index)+ " isn't contained in the structure.");
        Beads.push_back(bead{&Atoms[Index], NBead, Atoms[Index].Mass});
    }
    // if there's beads file, then coordinate of each bead can be set from the beginning
    if (strlen(Input.get_beads_file().c_str()))
    {
        ifstream InFileBeads(Input.get_beads_file().c_str(), ios::in);
        const size_t BeadNum = BeadIndexes.size();
        for (index_t k = 0 ; k < NBead ; ++k)
        {
            for (index_t i = 0 ; i < BeadNum ; ++i) InFileBeads >> Beads[i].Coords[k].x >> Beads[i].Coords[k].y >> Beads[i].Coords[k].z;
            InFileBeads.ignore(150, '\n');
        }
    }
}

const index_t cpihmc::box::find_ge_index(const std::vector<index_t> &Vec, const index_t Num)
{
    index_t Left = -1;
    index_t Right = Vec.size();
    while (Right - Left > 1)
    {
        index_t Mid = Left + (Right - Left) / 2;
        if (Vec[Mid] < Num) Left = Mid;
        else Right = Mid;
    };
    return Right;
}

const index_t cpihmc::box::element2atom(const index_t ElementIndex, const index_t ElementAtomIndex) const
{
    if (ElementIndex < 0) throw invalid_argument("The element index shouldn't be negative.");
    else if (ElementIndex >= NType) throw invalid_argument("The element index is out of the element number.");
    if (ElementAtomIndex < 0) throw invalid_argument("The atom index of the element shouldn't be negative.");
    else if (ElementAtomIndex >= Elements[ElementIndex].Number) throw invalid_argument("The atom index of the element is out of the atom number of the element.");
    index_t AtomIndex = ElementAtomIndex;
    for (index_t i = 0 ; i < ElementIndex ; ++i) AtomIndex += Elements[i].Number;
    return AtomIndex;
}

void cpihmc::box::update_kinetic_energy()
{
    KinEng = 0.0;
    for (const auto &i : AtomIndexes) KinEng += Atoms[i].Vel * Atoms[i].Vel * 0.5 * Atoms[i].Mass * MassScal;
}

void cpihmc::box::insert_atom(const index_t AtomIndex)
{
    if (AtomIndex < 0) throw invalid_argument("The atom index shouldn't be negative.");
    else if (AtomIndex >= NAtoms) throw invalid_argument("The atom index is out of the atom number of the structure.");
    vector<index_t>::const_iterator IterMask = find(Masks.begin(), Masks.end(), AtomIndex);
    if (IterMask != Masks.end()) Masks.erase(IterMask);
    else throw invalid_argument("The atom with the index "+to_string(AtomIndex)+" cannot be inserted.");
    const index_t GEIndex = find_ge_index(AtomIndexes, AtomIndex);
    AtomIndexes.insert(AtomIndexes.begin()+GEIndex, AtomIndex);
}

void cpihmc::box::insert_atom(const index_t ElementIndex, const index_t ElementAtomIndex){insert_atom(element2atom(ElementIndex, ElementAtomIndex));}

void cpihmc::box::remove_atom(const index_t AtomIndex)
{
    if (AtomIndex < 0) throw invalid_argument("The atom index shouldn't be negative.");
    else if (AtomIndex >= NAtoms) throw invalid_argument("The atom index is out of the atom number of the structure.");
    const index_t GEIndex = find_ge_index(Masks, AtomIndex);
    if (GEIndex < Masks.size() && Masks[GEIndex] == AtomIndex) throw invalid_argument("The atom with the index "+to_string(AtomIndex)+" cannot be removed.");
    else Masks.insert(Masks.begin()+GEIndex, AtomIndex);
    vector<index_t>::const_iterator IterAtom = find(AtomIndexes.begin(), AtomIndexes.end(), AtomIndex);
    AtomIndexes.erase(IterAtom);
    Atoms[AtomIndex].Vel.set_zero();
    Atoms[AtomIndex].Force.set_zero();
}

void cpihmc::box::remove_atom(const index_t ElementIndex, const index_t ElementAtomIndex){remove_atom(element2atom(ElementIndex, ElementAtomIndex));}