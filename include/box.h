//
// Created by Jin Bin on 2021/11/08.
//

#ifndef MD_MC_BOX_H
#define MD_MC_BOX_H

#include <vector>
#include "element.h"
#include "constant.h"
#include <random>
#include "matrix3.h"
#include <map>
using namespace std;

#define PI acos(-1)

extern constant Consts;

class box {
public:
    element *elements;
    map<string, int> element_type;
    Matrix3 lattice_vector;
    double *electron_number;
    int N_atoms;
    int ntype = Consts.ntype;
    double T = Consts.T;
    double kinetic_energy = 0;
    double quantum_kinetic_energy = 0;
    double potential_energy = 0;
    double energy = 0;
    const double k_B = 2.97116e-6;
    vector<bead> beads;
    vector<int> bead_index;
    box(element *, Matrix3, int, double *, long long);
    void init_velocities();
    void update_kinetic_energy();
    double Random() const;
    double *new_atoms(element *Elements) const;
    void reload_r(const double *vec) const;
    void regular() const; // move atom(s) outside the cell inside
    double RanGaussian(double mass) const;
private:
    long long seed;
    Matrix3 lattice_vector_inverse;
};


#endif //MD_MC_BOX_H
