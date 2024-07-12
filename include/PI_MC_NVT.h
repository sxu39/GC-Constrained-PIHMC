//
// Created by Jin Bin on 2022/02/21.
//

#ifndef MD_MC_PI_MC_NVT_H
#define MD_MC_PI_MC_NVT_H

#include "box.h"
#include <algorithm>
#include "NSphere.h"

class PI_MC_NVT {
public:
    box *Box;
    double (*calculate)(box *);
    int P = Consts.P;
    int J = Consts.J;
    int middle_atom = Consts.constraint_atom[1];
    int left_atom = Consts.constraint_atom[0];
    int after_atom = Consts.constraint_atom[2];
    atom *hidari_atom, *naka_atom, *migi_atom;
    double radial_step = Consts.radial;
    double theta_step = Consts.theta_step;
    double phi_step = Consts.phi_step;
    RanNSphere trial_move{radial_step};
    double react_coord = Consts.reaction_coordinate;
    explicit PI_MC_NVT(box *Box, double (*calculate)(box *));
    ~PI_MC_NVT();
    void evolve(int &current, bool centroid=false) const;
private:
    int *number_list;
    template <class T>
    static T min(T x, T y){
        return x<y?x:y;
    }
    void internal_evolve(int &current) const;
    void centroid_evolve(int &current) const;
    atom *change_atom_pointer() const;
    void u2x(bead *bead, int change_bead) const;
    // move the mass of center of beads to centroid after the change of bead or the change of centroid in two kinds of movement
    void back_centroid(bead *bead) const;
    void set_constraint_atoms();
    void scale(double &factor) const;
    void scale(atom *correlation_atom, atom *constraint_atom, double length) const;
};


#endif //MD_MC_PI_MC_NVT_H
