//
// Created by Jin Bin on 2022/02/21.
//

#ifndef CPIHMC_PI_HMC_NVT_H
#define CPIHMC_PI_HMC_NVT_H

#include "box.h"
#include "NSphere.h"
#include "wall.h"
#include <algorithm>

class PI_HMC_NVT {
public:
    box *Box;
    double (*calculate)(box *);
    int n_step = Consts.N_step; // the number of steps between two judge
    int P = Consts.P;
    int J = Consts.J;
    double delta_t = Consts.Delta_t;
    double mass_scale;
    int middle_atom = Consts.constraint_atom[1];
    int left_atom = Consts.constraint_atom[0];
    int after_atom = Consts.constraint_atom[2];
    atom *hidari_atom, *naka_atom, *migi_atom; // constraint atoms for difference situation
    atom *correlation_atom, *constraint_atom; // constraint atoms for distance situation
    double length_step = Consts.length_step;
    double theta_step = Consts.theta_step;
    double phi_step = Consts.phi_step;
    double react_coord;
    PI_HMC_NVT() = default;
    PI_HMC_NVT(box *Box, double (*calculate)(box *));
    virtual ~PI_HMC_NVT();
    void evolve(int &current, bool centroid=false, bool MC_step=false) const;
private:
    void u2x(bead *bead, int change_bead) const;
    double *adjust_angle(atom *move_atom, int type) const;
protected:
    hard_boundary *hardBoundary;
    // mass scaling coefficient used to scale the mass back for internal evolution
    const double PI_mass_scaling = 1 / Consts.M_scaling;
    template <class T>
    static T min(T x, T y){
        return x<y?x:y;
    }
    // move the mass of center of beads to centroid after the change of bead or the change of centroid in two kinds of movement
    void back_centroid(bead *bead) const;
    virtual void set_constraint_atoms();
    virtual void MC_evolve(int &current) const;
    virtual void HMC_evolve(int &current) const;
    void internal_evolve(int &current) const;
    void centroid_evolve(int &current, bool MC_step=false) const;
};

class PI_HMC_NVT_DIS : public PI_HMC_NVT{
public:
    PI_HMC_NVT_DIS() = default;
    PI_HMC_NVT_DIS(box *Box, double (*calculate)(box *));
protected:
    int correlation_index = Consts.constraint_atom[0]; // one atom of the constraint used as rigid body
    int constraint_index = Consts.constraint_atom[1]; // another atom of the constraint
    const double small_ball_r = Consts.radial;
    RanNSphere sphere{small_ball_r};
    virtual void MC_evolve(int &current) const;
    virtual void HMC_evolve(int &current) const;
    void set_constraint_atoms();
    static double calculate_sin_theta(Vector3<double> &vec);
    double adjust_angle() const;
};

#endif // CPIHMC_HMC_NVT_H
