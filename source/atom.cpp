//
// Created by Jin Bin on 2021/11/08.
//

#include "atom.h"

atom::atom():type(normal){}

void atom::force_initialize(){
    f.x = 0;
    f.y = 0;
    f.z = 0;
}

atom::~atom() = default;

bead::bead(atom *Atom, int P, double mass):bead_atom(Atom), mass(mass), kinetic_energy(0){
    r.resize(P);
    for (int i = 0 ; i < P ; ++i)
        r[i] = bead_atom->r;
}

bead::~bead() = default;