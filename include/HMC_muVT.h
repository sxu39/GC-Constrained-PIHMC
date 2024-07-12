//
// Created by Jin Bin on 2024/01/08
//

#ifndef CPIHMC_HMC_muVT_H
#define CPIHMC_HMC_muVT_H

#include "HMC_NVT.h"
#include "ml_potential.h"

extern int ne_reject;
extern int ne_num;

class HMC_muVT : public HMC_NVT {
public:
    HMC_muVT(box *, double (*)(box *));
    void evolve(int &current) const override;
private:
    double mu = Consts.mu * potential::eV2energy;
    pair<double, double> Ne_range = Consts.Ne_range;
    double delta_Ne = Consts.delta_Ne;
    void electron_number_evolve(int &) const;
};

#endif