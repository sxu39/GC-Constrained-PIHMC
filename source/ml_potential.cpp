//
// Created by Jin Bin on 2022/02/10.
//

#include "ml_potential.h"

deepmd::DeepPot *dp; // input the DP model
deepmd::DeepPot canonical_dp;
deepmd::DeepPot grand_canonical_dp;

double deep_potential::calculate(box *Box) {
    double potential_energy = 0;

    // obtain the index for bead atom
    vector<int> is_bead(Box->N_atoms, 0);
    for (int i = 0 ; i < Box->bead_index.size() ; ++i){
        is_bead[Box->bead_index[i]] = i+1;
        Box->beads[i].kinetic_energy = 0;
    }

    // collect all beads structures in the format required for DP inference
    vector<ValueType> coords(Consts.P*Box->N_atoms*3);
    vector<ValueType> cells;
    vector<int> atype(Box->N_atoms);
    provide_structures(coords, atype, cells, Box);

    // perform DP inference on all beads structures
    vector<double> energies;
    vector<ValueType> forces;

    if (strcmp(Consts.GC_type.c_str(), "Continuous") == 0){
        vector<ValueType> electron_numbers(Consts.P, *Box->electron_number);
        single_point_calculate(energies, forces, coords, atype, cells, electron_numbers);
    }
    else if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0)
        single_point_calculate(energies, forces, coords, atype, cells);

    // change position of each bead atom into certain bead situation and run cooresponding single point calculation and update force and energy
    vector<Vector3<double> > mean_force(Box->N_atoms);
    for (int i = 0 ; i < Consts.P ; ++i){
        potential_energy += eV2energy * energies[i] / Consts.P;
        for (int j = 0 ; j < Box->N_atoms ; ++j){
                Vector3<double> each_force;
                each_force.x = forces[(i*Box->N_atoms+j)*3];
                each_force.y = forces[(i*Box->N_atoms+j)*3+1];
                each_force.z = forces[(i*Box->N_atoms+j)*3+2];
                each_force *= eV_d_A2force / Consts.P;
                mean_force[j] += each_force;
                if (is_bead[j])
                    Box->beads[is_bead[j]-1].kinetic_energy += each_force * 0.5 * (Box->beads[is_bead[j]-1].bead_atom->r-Box->beads[is_bead[j]-1].r[i]);
        }
    }

    // update force for all atoms
    int curr = 0;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            Box->elements[i].atoms[j].f = mean_force[curr];
            ++curr;
        }
    }

    // update pseudo quantum kinetic energy
    Box->quantum_kinetic_energy = 0;
    for (auto &bead : Box->beads)
        Box->quantum_kinetic_energy += bead.kinetic_energy;
    return potential_energy;
}

void deep_potential::provide_structures(vector<ValueType> &coords, vector<int> &atype, vector<ValueType> &cells, box *Box){
    vector<ValueType > cell = {Box->lattice_vector.e11*Bohr2A, Box->lattice_vector.e12*Bohr2A, Box->lattice_vector.e13*Bohr2A,
                                Box->lattice_vector.e21*Bohr2A, Box->lattice_vector.e22*Bohr2A, Box->lattice_vector.e23*Bohr2A,
                                Box->lattice_vector.e31*Bohr2A, Box->lattice_vector.e32*Bohr2A, Box->lattice_vector.e33*Bohr2A};
    int current = 0;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j)
            atype[current+j] = Box->element_type[Box->elements[i].label];
        current += Box->elements[i].number;
    }

    vector<Vector3<double> > cent_pos(Box->beads.size()); // centroid position for bead atoms
    vector<Vector3<double> >::iterator cent_pos_iter;
    
    if (Consts.P > 1){
        // save centroid position for bead atoms
        cent_pos_iter = cent_pos.begin();
        for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
            *cent_pos_iter = it->bead_atom->r;
    }

    current = 0;
    for (int i = 0 ; i < Consts.P ; ++i){
        for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++)
            it->bead_atom->r = it->r[i];
        for (int j = 0 ; j < Box->ntype ; ++j){
            for (int k = 0 ; k < Box->elements[j].number ; ++k){
                coords[(current+k)*3] = Box->elements[j].atoms[k].r.x*Bohr2A;
                coords[(current+k)*3+1] = Box->elements[j].atoms[k].r.y*Bohr2A;
                coords[(current+k)*3+2] = Box->elements[j].atoms[k].r.z*Bohr2A;
            }
            current += Box->elements[j].number;
        }
        cells.insert(cells.end(), cell.begin(), cell.end());
    }

    if (Consts.P > 1){
        // recover the centroid position for bead atoms
        cent_pos_iter = cent_pos.begin();
        for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
            it->bead_atom->r = *cent_pos_iter;
    }
}

void deep_potential::single_point_calculate(vector<double> &energies, vector<ValueType> &forces, const vector<ValueType> &coords, const vector<int> &atype, const vector<ValueType> &cells) {
    vector<ValueType > v;
    dp->compute(energies, forces, v, coords, atype, cells);
}


void deep_potential::single_point_calculate(vector<double> &energies, vector<ValueType> &forces, const vector<ValueType> &coords, const vector<int> &atype, const vector<ValueType> &cells, const vector<ValueType> &electron_numbers) {
    vector<ValueType > v;
    dp->compute(energies, forces, v, coords, atype, cells, electron_numbers);
}
