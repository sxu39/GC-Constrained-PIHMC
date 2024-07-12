//
// Created by Jin Bin on 2022/02/21.
//

#ifndef MD_MC_HMC_NVT_H
#define MD_MC_HMC_NVT_H

#include "box.h"
#include "wall.h"

class HMC_NVT {
public:
    box *Box;
    double (*calculate)(box *);
    explicit HMC_NVT(box *Box, double (*calculate)(box *));
    ~HMC_NVT(){delete hardBoundary;}
    virtual void evolve(int &current) const;
private:
    hard_boundary *hardBoundary;
    int n_step = 3;
    double delta_t = Consts.Delta_t;
};


#endif //MD_MC_HMC_NVT_H
