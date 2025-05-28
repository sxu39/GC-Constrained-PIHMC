#include "atom.h"
using namespace cpihmc;

void cpihmc::bead::back_centroid()
{
    vec3<prec_t> displacement = BeadAtom->Coord;
    prec_t P_inv = 1/(prec_t)Coords.size();
    for (const auto &Coord : Coords)
        displacement -= Coord * P_inv;
    for (auto &Coord : Coords)
        Coord += displacement;
}