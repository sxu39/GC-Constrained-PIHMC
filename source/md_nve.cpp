#include "md_nve.h"
using namespace cpihmc;

void cpihmc::md_nve::evolve() const
{
    for (const auto &i : Box->AtomIndexes)
        if (Box->Atoms[i].Move.x * Box->Atoms[i].Move.y * Box->Atoms[i].Move.z)
        {
            Box->Atoms[i].Vel += Box->Atoms[i].Force * TimeStep * 0.5 / Box->Atoms[i].Mass;
            Box->Atoms[i].Coord += Box->Atoms[i].Vel * TimeStep;
        }
    Box->PotEng = Pot->infer(Box);
    for (const auto &i : Box->AtomIndexes)
        if (Box->Atoms[i].Move.x * Box->Atoms[i].Move.y * Box->Atoms[i].Move.z)
            Box->Atoms[i].Vel += Box->Atoms[i].Force * TimeStep * 0.5 / Box->Atoms[i].Mass;
    Box->update_total_energy(true);
}
