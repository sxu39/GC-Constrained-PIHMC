//
// Created by Jin Bin on 2022/08/18.
//

#ifndef CPIHMC_GEN_FUNC_H
#define CPIHMC_GEN_FUNC_H

#include "HMC_NVT.h"
#include "MC_NVT_RC.h"
#include "HMC_NVT_RC.h"
#include "PI_MC_NVT.h"
#include "PI_HMC_NVT.h"
#include "PI_HMC_muVT.h"
#include "cell.h"
#include <unordered_map>

double react_coord(atom *left, atom *right); // calculate reaction coordinate for distance situation
double react_coord(atom *left, atom *middle, atom *right); // calculate reaction coordinate
double react_coord(atom *left, atom *middle, atom *right, Matrix3 lattice); // calculate reaction coordinate considering mirror image in PBCs
void set_beads(box *Box); // set bead atom for box

void init_beads(box *Box, ifstream &ifbeads); // initialize the coordinates of beads according to file

void save_stru(int curr_step, cell *Cell, box *Box); // to save all structures of each iterative step
template <class T>
double internal_energy_estimator(T *mc) // calculate internal energy estimator
/*
    mc: MC simulation object
*/
{
    double internal_energy = 1.5 * mc->Box->N_atoms * mc->Box->k_B * mc->Box->T + mc->Box->quantum_kinetic_energy + mc->Box->potential_energy;
    return internal_energy;
}
template <class T>
double free_energy_estimator_left(T *mc); // calculate free energy estimator of left atom
template <class T>
double free_energy_estimator_right(T *mc) // calculate free energy estimator of right atom
/*
    mc: MC simulation object
*/
{
    double free_energy = 0;

    // derivative of potential energy term    
    Vector3<double> direction = mc->migi_atom->r - mc->naka_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy += direction * mc->migi_atom->f;

    // Jacobian term
    free_energy += 2*mc->Box->k_B*mc->Box->T * length_inv;
    return free_energy;
}
template <class T>
double free_energy_estimator_constraint(T *mc); // calculate free energy estimator of constraint atom
template <class T>
double free_energy_estimator_correlation(T *mc); // calculate free energy estimator of correlation atom
inline void update_physical_quantities(HMC_NVT *hmc, box *Box, unordered_map<string, double> &output_quantities, int i, double &ene_ave)
/* 
    update different physical quantities for each step in HMC method
    hmc: HMC simulation object
    Box: the simulated system
    output_quantities: mapping between output physical quantities and their label
    i: curr_step
    ene_ave: potential energy average
*/
{
    output_quantities["ke"] = Box->kinetic_energy;
    output_quantities["pe"] = Box->potential_energy;
    output_quantities["etotal"] = Box->energy;
    output_quantities["ne"] = *Box->electron_number;
    ene_ave = i / double(i+1) * ene_ave + Box->potential_energy / (i+1);
    output_quantities["pe_ave"] = ene_ave;
}
template <class T>
inline void update_physical_quantities(T *mc, box *Box, unordered_map<string, double> &output_quantities, int i, double &ene_ave, 
double &ene_2_ave)
/* 
    update different physical quantities for each step
    mc: MC simulation object
    Box: the simulated system
    output_quantities: mapping between output physical quantities and their label
    i: curr_step
    ene_ave: potential energy average
    ene_2_ave: the average of the square of potential energy
*/
{
    output_quantities["ke"] = Box->kinetic_energy;
    output_quantities["pe"] = Box->potential_energy;
    output_quantities["etotal"] = Box->energy;
    output_quantities["ne"] = *Box->electron_number;
    output_quantities["rc"] = react_coord(mc->hidari_atom, mc->naka_atom, mc->migi_atom);
    output_quantities["eint"] = internal_energy_estimator(mc);
    output_quantities["mfl"] = free_energy_estimator_left(mc);
    output_quantities["mfr"] = free_energy_estimator_right(mc);
    output_quantities["mf"] = 0.5 * (output_quantities["mfl"]+output_quantities["mfr"]);
    ene_ave = i / double(i+1) * ene_ave + Box->potential_energy / (i+1);
    ene_2_ave = i / double(i+1) * ene_2_ave + pow(Box->potential_energy, 2) / (i+1);
    output_quantities["pe_ave"] = ene_ave;
    output_quantities["pe_var"] = ene_2_ave - pow(ene_ave, 2);
}
template <class T>
inline void update_physical_quantities_distance(T *mc, box *Box, unordered_map<string, double> &output_quantities, int i, double &ene_ave, 
double &ene_2_ave)
/* 
    update different physical quantities for each step in distance situation
    mc: MC simulation object
    Box: the simulated system
    output_quantities: mapping between output physical quantities and their label
    i: curr_step
    ene_ave: potential energy average
    ene_2_ave: the average of the square of potential energy
*/
{
    output_quantities["ke"] = Box->kinetic_energy;
    output_quantities["pe"] = Box->potential_energy;
    output_quantities["etotal"] = Box->energy;
    output_quantities["ne"] = *Box->electron_number;
    output_quantities["rc"] = react_coord(mc->correlation_atom, mc->constraint_atom);
    output_quantities["eint"] = internal_energy_estimator(mc);
    output_quantities["mfl"] = free_energy_estimator_constraint(mc);
    output_quantities["mfr"] = free_energy_estimator_correlation(mc);
    output_quantities["mf"] = 0.5 * (output_quantities["mfl"]+output_quantities["mfr"]);
    ene_ave = i / double(i+1) * ene_ave + Box->potential_energy / (i+1);
    ene_2_ave = i / double(i+1) * ene_2_ave + pow(Box->potential_energy, 2) / (i+1);
    output_quantities["pe_ave"] = ene_ave;
    output_quantities["pe_var"] = ene_2_ave - pow(ene_ave, 2);
}
inline void output_physical_quantities(ofstream &output, unordered_map<string, int> &quantity_digits, unordered_map<string, double> &output_quantities)
/*
    output physical quantities required by the parameter file
    output: the out file stream
    quantity_digits: mapping between physical quantities' label and digits after point in the output number
    output_quantities: mapping between output physical quantities and their label
*/
{
    vector<string>::iterator it;
    for (it = Consts.output_terms.begin() ; it < Consts.output_terms.end()-1 ; ++it)
        output << setprecision(quantity_digits[*it]) << output_quantities[*it] << " ";
    output << setprecision(quantity_digits[*it]) << output_quantities[*it] << endl;
}
#endif