//
// Created by Jin Bin on 2022/02/10.
//

#include "first_principle.h"

extern cell Cell;

double vasp::calculate(box *Box) {
    double potential_energy = 0;
    vector<Vector3<double> > cent_pos(Box->beads.size()); // centroid position for bead atoms
    
    // save centroid position for bead atoms
    vector<Vector3<double> >::iterator cent_pos_iter;
    cent_pos_iter = cent_pos.begin();
    for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
        *cent_pos_iter = it->bead_atom->r;

    // folders to run sp calculation for each bead
    vector<string> folders(Consts.P);
    for (int i = 0 ; i < Consts.P ; ++i){
        char index[36];
        sprintf(index, format.c_str(), i);
        folders[i] = Consts.work_path + "_" + index;
    }
    
    // change position of each bead atom into certain bead situation and run cooresponding first principle calculation
    for (int i = 0 ; i < Consts.P ; ++i){
        for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++)
            it->bead_atom->r = it->r[i];
        string bead_atom_file;
        bead_atom_file = folders[i] + "/" + Consts.atom_file;
        Cell.print_stru_file(bead_atom_file);
        string command;
        // if (i == P-1)
        command = "cd " + folders[i] + " && stru2pos.py && mpirun -np 8 vasp_gam > log 2>err cd ..";
        // else
        //     command = "cd " + folders[i] + " && ABACUS > log 2>err & cd ..";
        system(command.c_str());
    }
    // system("sleep 6");

    // obtain energy and force from output, and calculate quantum kinetic energy via virial theorem
    vector<int> is_bead(Box->N_atoms, 0);
    for (int i = 0 ; i < Box->bead_index.size() ; ++i)
        is_bead[Box->bead_index[i]] = i+1;
    vector<Vector3<double> > mean_force(Box->N_atoms);
    for (int i = 0 ; i < Consts.P ; ++i){
        string command = "cd " + folders[i];
        command += " && ene > energy.dat";
        command += " && force > force.dat";
        command += " && cd ..";
        system(command.c_str());

        string ene_file = folders[i] + "/energy.dat";
        ifstream ife(ene_file.c_str(), ios::in);
        double temp_potential;
        if (ife.good())
            ife >> temp_potential;
        potential_energy += eV2energy * temp_potential / Consts.P;

        string force_file = folders[i] + "/force.dat";
        ifstream iff(force_file.c_str(), ios::in);
        if (iff.good()){
            for (int j = 0 ; j < Box->N_atoms ; ++j){
                Vector3<double> force;
                iff >> force.x >> force.y >> force.z;
                force *= eV_d_A2force / Consts.P;
                mean_force[j] += force;
                if (is_bead[j])
                    Box->beads[is_bead[j]-1].kinetic_energy += force * Box->beads[is_bead[j]-1].r[i] * 0.5;
                iff.ignore(150, '\n');
            }
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

    // recover the centroid position for bead atoms
    cent_pos_iter = cent_pos.begin();
    for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
        it->bead_atom->r = *cent_pos_iter;
    return potential_energy;
}

double abacus::calculate(box *Box) {
    double potential_energy = 0;
    vector<Vector3<double> > cent_pos(Box->beads.size()); // centroid position for bead atoms
    
    // save centroid position for bead atoms
    vector<Vector3<double> >::iterator cent_pos_iter;
    cent_pos_iter = cent_pos.begin();
    for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
        *cent_pos_iter = it->bead_atom->r;

    // folders to run sp calculation for each bead
    vector<string> folders(Consts.P);
    for (int i = 0 ; i < Consts.P ; ++i){
        char index[36];
        sprintf(index, format.c_str(), i);
        folders[i] = Consts.work_path + "_" + index;
    }
    
    // change position of each bead atom into certain bead situation and run cooresponding first principle calculation
    for (int i = 0 ; i < Consts.P ; ++i){
        for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++)
            it->bead_atom->r = it->r[i];
        string bead_atom_file;
        bead_atom_file = folders[i] + "/" + Consts.atom_file;
        Cell.print_stru_file(bead_atom_file);
        string command;
        // if (i == P-1)
        command = "cd " + folders[i] + " && ABACUS > log 2>err cd ..";
        // else
        //     command = "cd " + folders[i] + " && ABACUS > log 2>err & cd ..";
        system(command.c_str());
    }
    // system("sleep 6");

    // obtain energy and force from output, and calculate quantum kinetic energy via virial theorem
    vector<int> is_bead(Box->N_atoms, 0);
    for (int i = 0 ; i < Box->bead_index.size() ; ++i)
        is_bead[Box->bead_index[i]] = i+1;
    vector<Vector3<double> > mean_force(Box->N_atoms);
    for (int i = 0 ; i < Consts.P ; ++i){
        string command = "cd " + folders[i];
        command += " && extract_energy.sh > energy.dat";
        command += " && extract_force.sh > force.dat";
        command += " && cd ..";
        system(command.c_str());

        string ene_file = folders[i] + "/energy.dat";
        ifstream ife(ene_file.c_str(), ios::in);
        double temp_potential;
        if (ife.good())
            ife >> temp_potential;
        potential_energy += eV2energy * temp_potential / Consts.P;

        string force_file = folders[i] + "/force.dat";
        ifstream iff(force_file.c_str(), ios::in);
        if (iff.good()){
            for (int j = 0 ; j < Box->N_atoms ; ++j){
                Vector3<double> force;
                iff >> force.x >> force.y >> force.z;
                force *= eV_d_A2force / Consts.P;
                mean_force[j] += force;
                if (is_bead[j])
                    Box->beads[is_bead[j]-1].kinetic_energy += force * Box->beads[is_bead[j]-1].r[i] * 0.5;
                iff.ignore(150, '\n');
            }
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

    // recover the centroid position for bead atoms
    cent_pos_iter = cent_pos.begin();
    for (vector<bead>::iterator it = Box->beads.begin() ; it < Box->beads.end() ; it++, cent_pos_iter++)
        it->bead_atom->r = *cent_pos_iter;
    return potential_energy;
}