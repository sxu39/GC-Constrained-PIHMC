#include "gen_func.h"

double react_coord(atom *left, atom *right)
/*
    left: correlation atom of the constraint
    right: constraint atom of the constraint
*/
{
    Vector3<double> dis = right->r - left->r;
    return dis.norm();
}

double react_coord(atom *left, atom *middle, atom *right)
/*
    left: left atom of the constraint
    middle: middle atom of the constraint
    right: right atom of the constraint
*/
{
    Vector3<double> dis_1 = middle->r - left->r;
    Vector3<double> dis_2 = middle->r - right->r;
    return dis_1.norm() - dis_2.norm();
}

double react_coord(atom *left, atom *middle, atom *right, Matrix3 lattice)
/*
    left: left atom of the constraint
    middle: middle atom of the constraint
    right: right atom of the constraint
    lattice: lattice vectors of the cell, only for zero off-diagonal elements
*/
{
    Vector3<double> dis_1 = middle->r - left->r;
    Vector3<double> dis_2 = middle->r - right->r;

    // adjust distance between two atoms in order to make each direction of it less than half of the lattice constant
    dis_1.x -= lattice.e11 * floor(dis_1.x/lattice.e11+0.5);
    dis_1.y -= lattice.e22 * floor(dis_1.y/lattice.e22+0.5);
    dis_1.z -= lattice.e33 * floor(dis_1.z/lattice.e33+0.5);
    dis_2.x -= lattice.e11 * floor(dis_2.x/lattice.e11+0.5);
    dis_2.y -= lattice.e22 * floor(dis_2.y/lattice.e22+0.5);
    dis_2.z -= lattice.e33 * floor(dis_2.z/lattice.e33+0.5);
    return dis_1.norm() - dis_2.norm();
}

void set_beads(box *Box)
/*
    Box: the simulated system
*/
{
    vector<int> bead_index = Consts.bead_index;
    sort(bead_index.begin(), bead_index.end());
    int current = 0;
    int index = 0;
    int size = bead_index.size();
    vector<pair<int, int>> bead_atom_index;
    Box->bead_index = bead_index;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            if (index == size)
                break;
            if (current == bead_index[index]){
                bead new_bead(&(Box->elements[i].atoms[j]), Consts.P, Box->elements[i].mass);
                Box->beads.push_back(new_bead);
                bead_atom_index.push_back(pair<int, int>(i, j));
                ++index;
            }
            ++current;
        }
        if (index == size)
            break;
    }
    for (index = 0 ; index < size ; ++index)
        Box->elements[bead_atom_index[index].first].atoms[bead_atom_index[index].second].pbead = &(Box->beads[index]);
}

void init_beads(box *Box, ifstream &ifbeads)
/*
    Box: the simulated system
    ifbeads: beads coordinates file input stream
*/
{
    int bead_num = Box->bead_index.size();
    for (int k = 0 ; k < Consts.P ; ++k){
        for (int i = 0 ; i < bead_num ; ++i)
            ifbeads >> Box->beads[i].r[k].x >> Box->beads[i].r[k].y >> Box->beads[i].r[k].z;
        ifbeads.ignore(150, '\n');
    }
}

void save_stru(int curr_step, cell *Cell, box *Box)
/*
    curr_step: step index
    cell: cell from structure file
    Box: box system
*/
{
    // output current structure with centroid structure file and beads file
    Cell->print_stru_file("ALL_STRU/" + Consts.atom_file + "_" + to_string(curr_step));
    if (Consts.P > 1)
    // output current beads for quantum situations
    {
        string beads_file = "ALL_STRU/" + Consts.beads_file + "_" + to_string(curr_step);
        ofstream ofbeads(beads_file.c_str());
        vector<bead>::iterator it;
        for (int i = 0 ; i < Consts.P ; ++i){
            for (it = Box->beads.begin() ; it < Box->beads.end() - 1; it++)
                ofbeads << it->r[i] << "\t";
            ofbeads << it->r[i] << endl;
        }
        ofbeads.close();
    }
}

template double free_energy_estimator_left<MC_NVT_RC>(MC_NVT_RC *mc);
template double free_energy_estimator_left<HMC_NVT_RC>(HMC_NVT_RC *mc);
template double free_energy_estimator_left<PI_MC_NVT>(PI_MC_NVT *mc);
template double free_energy_estimator_left<PI_HMC_NVT>(PI_HMC_NVT *mc);
template double free_energy_estimator_left<PI_HMC_muVT>(PI_HMC_muVT *mc);

template <class T>
double free_energy_estimator_left(T *mc)
/*
    mc: MC simulation object
*/
{
    double free_energy = 0;

    // derivative of potential energy term
    Vector3<double> direction = mc->hidari_atom->r - mc->naka_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy -= direction * mc->hidari_atom->f;

    // Jacobian term
    free_energy -= 2*mc->Box->k_B*mc->Box->T * length_inv;
    return free_energy;
}

template double free_energy_estimator_constraint<MC_NVT_DIS>(MC_NVT_DIS *mc);
template double free_energy_estimator_constraint<HMC_NVT_DIS>(HMC_NVT_DIS *mc);
template double free_energy_estimator_constraint<HMC_NVT_RC>(HMC_NVT_RC *mc);
template double free_energy_estimator_constraint<PI_HMC_NVT_DIS>(PI_HMC_NVT_DIS *mc);
template double free_energy_estimator_constraint<PI_HMC_NVT>(PI_HMC_NVT *mc);
template double free_energy_estimator_constraint<PI_HMC_muVT>(PI_HMC_muVT *mc);

template <class T>
double free_energy_estimator_constraint(T *mc)
/*
    mc: MC simulation object
*/
{
    double free_energy = 0;

    // derivative of potential energy term
    Vector3<double> direction = mc->constraint_atom->r - mc->correlation_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy -= direction * mc->constraint_atom->f;

    // Jacobian term
    free_energy -= 2*mc->Box->k_B*mc->Box->T * length_inv;
    return free_energy;
}

template double free_energy_estimator_correlation<MC_NVT_DIS>(MC_NVT_DIS *mc);
template double free_energy_estimator_correlation<HMC_NVT_DIS>(HMC_NVT_DIS *mc);
template double free_energy_estimator_correlation<HMC_NVT_RC>(HMC_NVT_RC *mc);
template double free_energy_estimator_correlation<PI_HMC_NVT_DIS>(PI_HMC_NVT_DIS *mc);
template double free_energy_estimator_correlation<PI_HMC_NVT>(PI_HMC_NVT *mc);
template double free_energy_estimator_correlation<PI_HMC_muVT>(PI_HMC_muVT *mc);

template <class T>
double free_energy_estimator_correlation(T *mc)
/*
    mc: MC simulation object
*/
{
    double free_energy = 0;

    // derivative of potential energy term
    Vector3<double> direction = mc->correlation_atom->r - mc->constraint_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy -= direction * mc->correlation_atom->f;

    // Jacobian term
    free_energy -= 2*mc->Box->k_B*mc->Box->T * length_inv;
    return free_energy;
}