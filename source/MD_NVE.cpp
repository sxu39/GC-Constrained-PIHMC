//
// Created by Jin Bin on 2022/02/21.
//

#include "MD_NVE.h"

MD_NVE::MD_NVE(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    Box->energy = Box->kinetic_energy + Box->potential_energy;
}

void MD_NVE::evolve(double delta_t) const{
    for (int j = 0 ; j < Box->ntype ; ++j){
        element *curr_ele = &(Box->elements[j]);
        for (int k = 0 ; k < curr_ele->number ; ++k){
            if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t / 2) / curr_ele->mass;
                curr_ele->atoms[k].r += curr_ele->atoms[k].v * delta_t;
            }
        }
    }
    Box->potential_energy = calculate(Box);
    for (int j = 0 ; j < Box->ntype ; ++j){
        element *curr_ele = &(Box->elements[j]);
        for (int k = 0 ; k < curr_ele->number ; ++k)
            if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z)
                curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t / 2) / curr_ele->mass;
    }
    Box->regular();
    Box->update_kinetic_energy();
    Box->energy = Box->kinetic_energy + Box->potential_energy;
}
