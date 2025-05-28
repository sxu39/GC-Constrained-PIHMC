#include "dof.h"
using namespace cpihmc;
using namespace std;

const vector<vec3<prec_t> *> cpihmc::obtain_coord_pointer(const vector<atom *> &AllAtoms)
{
    vector<vec3<prec_t> *> Coords;
    for (const auto &Atom : AllAtoms) Coords.push_back(&Atom->Coord);
    return Coords;
}

const prec_t cpihmc::obtain_mean_mass(const vector<atom *> &AllAtoms)
{
    prec_t MassSum = 0.0;
    for (const auto &Atom : AllAtoms) MassSum += Atom->Mass;
    return MassSum / AllAtoms.size();
}

cpihmc::dof_atom::dof_atom(const vector<atom *> &AllAtoms):AllAtoms(AllAtoms), Coord(obtain_coord_pointer(AllAtoms)), Mass(obtain_mean_mass(AllAtoms)){}

void cpihmc::dof_atom::update_force()
{
    Force.set_zero();
    for (const auto &Atom : AllAtoms) Force += Atom->Force;
}

cpihmc::dof::dof(atom &Atom):VirtAtom(nullptr), Coord(Atom.Coord), Vel(Atom.Vel), Force(Atom.Force), Mass(Atom.Mass), AtomNum(1){}

cpihmc::dof::dof(dof_atom &VirtAtom):VirtAtom(&VirtAtom), Coord(VirtAtom.Coord), Vel(VirtAtom.Vel), Force(VirtAtom.Force), Mass(VirtAtom.Mass), AtomNum(VirtAtom.AllAtoms.size()){}

const vec3<prec_t> cpihmc::dof::get_force() const
{
    if (VirtAtom) VirtAtom->update_force();
    return Force;
}

void cpihmc::dof::synchronize() const {if (VirtAtom) VirtAtom->Coord.synchronize();}