//
// Created by Jin Bin on 2021/11/08.
//

#ifndef MD_MC_ATOM_H
#define MD_MC_ATOM_H

#include <cmath>
#include "vector3.h"
#include <vector>

enum constraint {normal, hidari, naka, migi};

class bead;
class atom {
public:
    Vector3<double> r;
    Vector3<double> v = Vector3<double>(0, 0, 0);
    Vector3<double> f;
    Vector3<bool> move;
    constraint type;
    bead *pbead = nullptr;
    atom();
    void force_initialize();
    ~atom();
};

class bead {
public:
    atom * bead_atom;
    vector<Vector3<double> > r;
    double kinetic_energy;
    double mass;
    bead(atom *Atom, int P, double mass);
    ~bead();
};


#endif //MD_MC_ATOM_H
