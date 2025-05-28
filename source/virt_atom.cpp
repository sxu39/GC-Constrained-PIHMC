#include "virt_atom.h"
using namespace cpihmc;
using namespace std;

cpihmc::virt_atom::virt_atom(const vector<atom *> &AllAtoms, const prec_t Temp, const prec_t TimeStep, const prec_t MassScal, rand &Rand):AtomNum(AllAtoms.size()), AllAtoms(AllAtoms), Centroid(AllAtoms), TFAtoms(AtomNum+1, atom{obtain_mean_mass(AllAtoms)}), hmc(Temp, TimeStep, MassScal, Rand)
{
    const prec_t Mass = AllAtoms[0]->Mass;
    for (index_t i = 0 ; i < AtomNum ; ++i) if (AllAtoms[i]->Mass != Mass) throw invalid_argument("The mass of the atoms for the DOF should be the same.");
    update_virt_atom_coordinates();
}

void cpihmc::virt_atom::update_virt_atom_coordinates(){for (index_t Index = 1 ; Index < AtomNum ; ++Index) TFAtoms[Index].Coord = AllAtoms[Index-1]->Coord - AllAtoms[Index]->Coord;}

const vec3<prec_t> cpihmc::virt_atom::get_transformed_force(const index_t Index) const
{
    vec3<prec_t> TFForce;
    for (index_t i = 0 ; i < AtomNum ; ++i) TFForce += AllAtoms[i]->Force * ((prec_t)AtomNum * (i<Index) - Index);
    return TFForce / (prec_t)AtomNum;
}

void cpihmc::virt_atom::random_velocities()
{
    vec3<prec_t> Vel;
    for (index_t Index = 1 ; Index < AtomNum ; ++Index)
    {
        Vel.set_zero();
        for (index_t i = 0 ; i < AtomNum ; ++i) Vel += AllAtoms[i]->Vel * ((prec_t)AtomNum * (i<Index) - Index);
        Vel /= (prec_t)AtomNum;
        TFAtoms[Index].Vel = Vel;
    }
}

void cpihmc::virt_atom::update_atom_coordinates() const
{
    vec3<prec_t> Coord;
    for (index_t i = 0 ; i < AtomNum ; ++i)
    {
        Coord = Centroid.Coord;
        for (index_t Index = 1 ; Index < AtomNum ; ++Index) Coord += TFAtoms[Index].Coord * ((prec_t)AtomNum * (i<Index) - Index) / (prec_t)AtomNum;
        AllAtoms[i]->Coord = Coord;
    }
}

void cpihmc::virt_atom::hmc_evolve_forth_before()
{
    for (index_t Index = 1 ; Index < AtomNum ; ++Index) TFAtoms[Index].Vel += get_transformed_force(Index) * TimeStep * 0.5 / (TFAtoms[Index].Mass * MassScal);
    for (index_t Index = 1 ; Index < AtomNum ; ++Index) TFAtoms[Index].Coord += (-TFAtoms[Index-1].Vel + 2.0 * TFAtoms[Index].Vel - TFAtoms[Index+1].Vel) * TimeStep;
    update_atom_coordinates();
}

void cpihmc::virt_atom::hmc_evolve_forth_after()
{for (index_t Index = 1 ; Index < AtomNum ; ++Index) TFAtoms[Index].Vel += get_transformed_force(Index) * TimeStep * 0.5 / (TFAtoms[Index].Mass * MassScal);}

void cpihmc::virt_atom::hmc_evolve_back(){update_virt_atom_coordinates();}

void cpihmc::virt_atom::update_atom_velocities() const
{for (index_t i = 0 ; i < AtomNum ; ++i) AllAtoms[i]->Vel = Centroid.Vel + TFAtoms[i+1].Vel - TFAtoms[i].Vel;}