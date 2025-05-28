#include "ml_pot.h"
using namespace cpihmc;
using namespace std;

prec_t cpihmc::deep_pot::infer(box *Box)
{
    prec_t PotEng = 0;

    // obtain the index for bead atom
    vector<index_t> IsBead(Box->NAtoms, 0);
    for (index_t i = 0 ; i < Box->BeadIndexes.size() ; ++i)
    {
        IsBead[Box->BeadIndexes[i]] = i+1;
        Box->Beads[i].KinEng = 0.0;
    }

    // collect all beads structures in the format required for DP inference
    const index_t RealAtomNum = Box->AtomIndexes.size();
    vector<prec_t> Coords(Box->NBead*RealAtomNum*3);
    vector<prec_t> Cells;
    vector<index_t> AtomTypes(RealAtomNum);
    provide_structures(Coords, AtomTypes, Cells, Box);

    // perform DP inference on all beads structures
    vector<prec_t> Energies;
    vector<prec_t> Forces;
    if (DeepPot.dim_fparam()) single_point_calculate(Energies, Forces, Coords, AtomTypes, Cells, vector<prec_t>(Box->NBead, Box->ElecNum));
    else single_point_calculate(Energies, Forces, Coords, AtomTypes, Cells);

    // change position of each bead atom into certain bead situation and run cooresponding single point calculation and update force and energy
    const vector<index_t> AtomIndexes = Box->AtomIndexes;
    vector<vec3<prec_t>> MeanForce(RealAtomNum);
    for (index_t i = 0 ; i < Box->NBead ; ++i)
    {
        PotEng += eV2Energy * Energies[i] / Box->NBead;
        for (index_t j = 0 ; j < RealAtomNum ; ++j)
        {
            vec3<prec_t> EachForce;
            EachForce.x = Forces[(i*RealAtomNum+j)*3];
            EachForce.y = Forces[(i*RealAtomNum+j)*3+1];
            EachForce.z = Forces[(i*RealAtomNum+j)*3+2];
            EachForce *= eV_d_Ang2Force / Box->NBead;
            MeanForce[j] += EachForce;
            if (IsBead[AtomIndexes[j]])
                Box->Beads[IsBead[AtomIndexes[j]]-1].KinEng += EachForce * 0.5 * (Box->Beads[IsBead[AtomIndexes[j]]-1].BeadAtom->Coord-Box->Beads[IsBead[AtomIndexes[j]]-1].Coords[i]);
        }
    }

    // update force for all atoms
    for (index_t i = 0 ; i < RealAtomNum ; ++i)
        Box->Atoms[AtomIndexes[i]].Force = MeanForce[i];

    // update pseudo quantum kinetic energy
    Box->QtmKinEng = 0.0;
    for (const auto &Bead : Box->Beads)
        Box->QtmKinEng += Bead.KinEng;
    return PotEng;
}

void cpihmc::deep_pot::provide_structures(vector<prec_t> &Coords, vector<index_t> &AtomTypes, vector<prec_t> &Cells, box *Box)
{
    vector<prec_t> Cell = {Box->LatVec.e11*Bohr2Ang, Box->LatVec.e12*Bohr2Ang, Box->LatVec.e13*Bohr2Ang,
                           Box->LatVec.e21*Bohr2Ang, Box->LatVec.e22*Bohr2Ang, Box->LatVec.e23*Bohr2Ang,
                           Box->LatVec.e31*Bohr2Ang, Box->LatVec.e32*Bohr2Ang, Box->LatVec.e33*Bohr2Ang};
    const vector<index_t> &AtomIndexes = Box->AtomIndexes;
    const size_t RealAtomNum = AtomIndexes.size();
    index_t ElementIndex = 0;
    index_t CurrNum = 0;
    index_t Type;
    for (index_t i = 0 ; i < RealAtomNum ; ++i)
    {
        if (AtomIndexes[i] >= CurrNum)
        {
            CurrNum += Box->Elements[ElementIndex].Number;
            Type = Box->ElementType.at(Box->Elements[ElementIndex].Label);
            ++ElementIndex;
        }
        AtomTypes[i] = Type;
    }

    vector<vec3<prec_t>> CentPos(Box->Beads.size()); // centroid position for bead atoms
    vector<vec3<prec_t>>::iterator CentPosIter;
    
    if (Box->NBead > 1)
    {
        // save centroid position for bead atoms
        CentPosIter = CentPos.begin();
        for (vector<bead>::iterator it = Box->Beads.begin() ; it < Box->Beads.end() ; ++it, ++CentPosIter)
            *CentPosIter = it->BeadAtom->Coord;
    }

    for (index_t i = 0 ; i < Box->NBead ; ++i)
    {
        for (vector<bead>::iterator it = Box->Beads.begin() ; it < Box->Beads.end() ; ++it)
            it->BeadAtom->Coord = it->Coords[i];
        for (index_t j = 0 ; j < RealAtomNum ; ++j)
        {
            Coords[3*(i*RealAtomNum+j)] = Box->Atoms[AtomIndexes[j]].Coord.x * Bohr2Ang;
            Coords[3*(i*RealAtomNum+j)+1] = Box->Atoms[AtomIndexes[j]].Coord.y * Bohr2Ang;
            Coords[3*(i*RealAtomNum+j)+2] = Box->Atoms[AtomIndexes[j]].Coord.z * Bohr2Ang;
        }
        Cells.insert(Cells.end(), Cell.begin(), Cell.end());
    }

    if (Box->NBead > 1)
    {
        // recover the centroid position for bead atoms
        CentPosIter = CentPos.begin();
        for (vector<bead>::iterator it = Box->Beads.begin() ; it < Box->Beads.end() ; ++it, ++CentPosIter)
            it->BeadAtom->Coord = *CentPosIter;
    }
}

void cpihmc::deep_pot::single_point_calculate(vector<prec_t> &Energies, vector<prec_t> &Forces, const vector<prec_t> &Coords, const vector<index_t> &AtomTypes, const vector<prec_t> &Cells)
{
    vector<prec_t> Virials;
    DeepPot.compute(Energies, Forces, Virials, Coords, AtomTypes, Cells);
}

void cpihmc::deep_pot::single_point_calculate(vector<prec_t> &Energies, vector<prec_t> &Forces, const vector<prec_t> &Coords, const vector<index_t> &AtomTypes, const vector<prec_t> &Cells, const vector<prec_t> &ElecNums)
{
    vector<prec_t> Virials;
    DeepPot.compute(Energies, Forces, Virials, Coords, AtomTypes, Cells, ElecNums);
}