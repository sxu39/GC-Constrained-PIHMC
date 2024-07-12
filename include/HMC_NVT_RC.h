//
// Created by Jin Bin on 2022/02/21.
//

#ifndef CPIHMC_HMC_NVT_RC_H
#define CPIHMC_HMC_NVT_RC_H

#include "box.h"
#include "NSphere.h"
#include "wall.h"

class HMC_NVT_RC {
public:
    box *Box;
    double (*calculate)(box *);
    int n_step = Consts.N_step;
    int middle_atom = Consts.constraint_atom[1];
    int left_atom = Consts.constraint_atom[0];
    int after_atom = Consts.constraint_atom[2];
    atom *hidari_atom, *naka_atom, *migi_atom;
    atom *correlation_atom, *constraint_atom;
    double delta_t = Consts.Delta_t;
    double react_coord;
    double mass_scale;
    double length_step = Consts.length_step;
    double theta_step = Consts.theta_step;
    double phi_step = Consts.phi_step;
    explicit HMC_NVT_RC(box *Box, double (*calculate)(box *));
    HMC_NVT_RC() = default;
    virtual ~HMC_NVT_RC(){
        delete hardBoundary;
    }
    virtual void evolve(int &current, bool MC_step=false) const; /* if this function is labeled by const, then member variables of the class can also
    be labeled as const, then I can only define const reference, such as : atom * const & curr_atom = hidari_atom */
private:
    void set_constraint_atoms();
    // set atoms related to constraint before the evolution, include defining pointers to these atoms, set their types and
    // calculation scale factor for mass
    double *adjust_angle(atom *move_atom, int type) const;
protected:
    hard_boundary *hardBoundary;
    template <class T>
    static T min(T x, T y){
        return x<y?x:y;
    }
    static double calculate_sin_theta(Vector3<double> vec);
};

class HMC_NVT_DIS : public HMC_NVT_RC {
public:
    HMC_NVT_DIS(box *Box, double (*calculate)(box *));
    void evolve(int &current, bool MC_step=false) const;
private:
    int correlation_index = Consts.constraint_atom[0]; // one atom of the constraint used as rigid body
    int constraint_index = Consts.constraint_atom[1]; // another atom of the constraint
    const double small_ball_r = Consts.radial;
    RanNSphere sphere{small_ball_r};
    void MC_evolve(int &current) const;
    void HMC_evolve(int &current) const;
    void set_constraint_atoms();
    double adjust_angle() const;
};

#endif //CPIHMC_HMC_NVT_RC_H
