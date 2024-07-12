//
// Created by Jin Bin on 2023/06/26.
//

#ifndef CPIHMC_PI_HMC_muVT_H
#define CPIHMC_PI_HMC_muVT_H

#include "PI_HMC_NVT.h"
#include "ml_potential.h"

extern int ne_reject;

class PI_HMC_muVT : public PI_HMC_NVT_DIS {
public:
    PI_HMC_muVT(box *, double (*)(box *));
    void evolve(int &, bool Ne=false, bool centroid=false, bool MC_step=false) const;
private:
    double mu = Consts.mu * potential::eV2energy;
    pair<double, double> Ne_range = Consts.Ne_range;
    double delta_Ne = Consts.delta_Ne;
    vector<deepmd::DeepPot> dp_models;
    void init_models();
    void MC_evolve(int &current) const;
    void HMC_evolve(int &current) const;
    void electron_number_evolve(int &) const;
    void set_model() const;
};

#endif