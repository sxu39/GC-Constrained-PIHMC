//
// Created by Jin Bin on 2021/11/03.
//

#ifndef MD_NVE_LJ_MC_NVT_RC_H
#define MD_NVE_LJ_MC_NVT_RC_H

#include "box.h"
#include "NSphere.h"
#include "matrix3.h"
// #include "constant.h"
// #include <vector>
// #include <algorithm>
// #include <random>

extern constant Consts;

class MC_NVT_RC {
public:
    box *Box;
    double (*calculate)(box *);
    double reaction_coordinate;
    double radial_step = Consts.radial;
    double theta_step = Consts.theta_step;
    double phi_step = Consts.phi_step;
    RanNSphere trial_move{radial_step};
    int middle_atom = Consts.constraint_atom[1];
    int left_atom = Consts.constraint_atom[0];
    int after_atom = Consts.constraint_atom[2];
    atom *hidari_atom, *naka_atom, *migi_atom;
    atom *correlation_atom, *constraint_atom;
    explicit MC_NVT_RC(box *Box, double (*calculate)(box *));
    MC_NVT_RC() = default;
    void evolve(int &current) const;
    double free_energy_estimator_left() const;
    double free_energy_estimator_right() const;
    virtual ~MC_NVT_RC();
private:
    void set_constraint_atoms();
    void scale(double &factor) const;
protected:
    template <class T>
    static T min(T x, T y);
    int *number_list;
    atom *change_atom_pointer() const;
    void scale(atom *correlation_atom, atom *constraint_atom, double length) const;
};

class MC_NVT_DIS : public MC_NVT_RC {
public:
    MC_NVT_DIS(box *Box, double (*calculate)(box *));
    void evolve(int &current) const;
private:
    int correlation_index = Consts.constraint_atom[0]; // one atom of the constraint used as rigid body
    int constraint_index = Consts.constraint_atom[1]; // another atom of the constraint
    void set_constraint_atoms();
};

#endif //MD_NVE_LJ_MC_NVT_H
