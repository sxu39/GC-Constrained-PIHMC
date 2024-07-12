//
// Created by Jin Bin on 2022/02/21.
//

#include "HMC_NVT.h"

HMC_NVT::HMC_NVT(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    hardBoundary = new wall {Box};
}

void HMC_NVT::evolve(int &current) const {
    auto *temp_atoms = Box->new_atoms(Box->elements);
    double temp_energy = Box->energy;
    for (int step = 0 ; step < n_step; ++step) {
        for (int j = 0 ; j < Box->ntype ; ++j){
            element *curr_ele = &(Box->elements[j]);
            for (int k = 0 ; k < curr_ele->number ; ++k)
                if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                    curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t * 0.5) / curr_ele->mass;
                    curr_ele->atoms[k].r += curr_ele->atoms[k].v * delta_t;
                }
        }
        Box->potential_energy = calculate(Box);
        for (int j = 0 ; j < Box->ntype ; ++j){
            element *curr_ele = &(Box->elements[j]);
            for (int k = 0 ; k < curr_ele->number ; ++k)
                if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z)
                    curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t * 0.5) / curr_ele->mass;
        }
    }
    Box->update_kinetic_energy();
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    if (Box->Random() > min(1.0, exp(-(Box->energy-temp_energy)/Box->k_B/Box->T)) || !(hardBoundary->judge_accept())){
        Box->reload_r(temp_atoms);
        Box->potential_energy = calculate(Box);
        current++;
    }
    Box->init_velocities();
    Box->update_kinetic_energy();
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    delete [] temp_atoms;
}
