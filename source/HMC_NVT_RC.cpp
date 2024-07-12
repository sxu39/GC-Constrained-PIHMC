//
// Created by Jin Bin on 2022/02/21.
//

#include "HMC_NVT_RC.h"
#include "gen_func.h"

enum MC_evolve {length, theta_1, phi_1, theta_2, phi_2}; // the type of MC step

HMC_NVT_RC::HMC_NVT_RC(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    set_constraint_atoms();
    if (!Consts.set_rc){
        Consts.reaction_coordinate = ::react_coord(hidari_atom, naka_atom, migi_atom);
    }
    react_coord = Consts.reaction_coordinate;
    hardBoundary = new wall {Box};
}

void HMC_NVT_RC::evolve(int &current, bool MC_step) const {
    if (MC_step){
        double temp_energy = Box->potential_energy; // record current potential energy
        int kind = floor(Box->Random()*5); // choose MC step type at random
        if (kind == length){
            Vector3<double> left_dis = hidari_atom->r - naka_atom->r;
            Vector3<double> right_dis = migi_atom->r - naka_atom->r;
            double curr_dis = right_dis.norm();
            double new_length = curr_dis + (Box->Random()-0.5) * length_step;
            hidari_atom->r = naka_atom->r + (new_length+react_coord) * left_dis / left_dis.norm();
            migi_atom->r = naka_atom->r + new_length * right_dis / right_dis.norm();
            if (hardBoundary->judge_accept()){    
                Box->potential_energy = calculate(Box);
                double factor = pow((new_length+react_coord)*new_length/curr_dis/(curr_dis+react_coord), 2);
                if (Box->Random() > min(1.0, factor*exp(-(Box->potential_energy-temp_energy)/(Box->k_B*Box->T)))){
                    hidari_atom->r = naka_atom->r + (curr_dis+react_coord) * left_dis / left_dis.norm();
                    migi_atom->r = naka_atom->r + curr_dis * right_dis / right_dis.norm();
                    Box->potential_energy = calculate(Box);
                    current++;
                }
            }
            else {
                hidari_atom->r = naka_atom->r + (curr_dis+react_coord) * left_dis / left_dis.norm();
                migi_atom->r = naka_atom->r + curr_dis * right_dis / right_dis.norm();
                current++;
            }
            Box->energy = Box->kinetic_energy + Box->potential_energy;
        }
        else{
            atom *curr_atom = hidari_atom;
            if (kind > 2)
                curr_atom = migi_atom;
            Vector3<double> temp_pos = curr_atom->r;
            double *sin_theta = adjust_angle(curr_atom, kind);
            if (hardBoundary->judge_accept()){
                Box->potential_energy = calculate(Box);
                double factor = 1;
                if (kind % 2)
                    factor = sin_theta[1]/sin_theta[0];
                if (Box->Random() > min(1.0, factor*exp(-(Box->potential_energy-temp_energy)/(Box->k_B*Box->T)))){
                    curr_atom->r = temp_pos;
                    Box->potential_energy = calculate(Box);
                    current++;
                }
            }
            else {
                curr_atom->r = temp_pos;
                current++;
            }
            Box->energy = Box->kinetic_energy + Box->potential_energy;
            delete [] sin_theta;
        } 
    }
    else{
        double temp_energy = Box->energy;
        auto *temp_atoms = Box->new_atoms(Box->elements);
        naka_atom->v *= (1/pow(mass_scale, 0.5)); // velocity is sampled for single atom, so it's necessary to scale for the rigid part.
        for (int step = 0 ; step < n_step; ++step) {
            for (int j = 0 ; j < Box->ntype ; ++j){
                element *curr_ele = &(Box->elements[j]);
                for (int k = 0 ; k < curr_ele->number ; ++k){
                    if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                        if (curr_ele->atoms[k].type == naka)
                            curr_ele->atoms[k].v += (naka_atom->f+hidari_atom->f+migi_atom->f) * (delta_t*0.5) / (curr_ele->mass*mass_scale);
                        else if (curr_ele->atoms[k].type == normal)
                            curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t*0.5) / curr_ele->mass;
                    }
                }
            }
            for (int j = 0 ; j < Box->ntype ; ++j){
                element *curr_ele = &(Box->elements[j]);
                for (int k = 0 ; k < curr_ele->number ; ++k){
                    if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                        if (curr_ele->atoms[k].type == hidari || curr_ele->atoms[k].type == migi)
                            curr_ele->atoms[k].r += naka_atom->v * delta_t;
                        else
                            curr_ele->atoms[k].r += curr_ele->atoms[k].v * delta_t;
                    }
                }
            }
            Box->potential_energy = calculate(Box);
            for (int j = 0 ; j < Box->ntype ; ++j){
                element *curr_ele = &(Box->elements[j]);
                for (int k = 0 ; k < curr_ele->number ; ++k)
                    if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                        if (curr_ele->atoms[k].type == naka)
                            curr_ele->atoms[k].v += (naka_atom->f+hidari_atom->f+migi_atom->f) * (delta_t*0.5) / (curr_ele->mass*mass_scale);
                        else if (curr_ele->atoms[k].type == normal)
                            curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t*0.5) / curr_ele->mass;
                    }
            }
        }
        naka_atom->v *= pow(mass_scale, 0.5); // scale the velocity back to calculation kinetic energy
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
}

void HMC_NVT_RC::set_constraint_atoms(){
    int current = 0;
    int max_index = max(left_atom, max(middle_atom, after_atom));
    mass_scale = 0;
    double naka_mass;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            if (current > max_index)
                break;
            if (current == left_atom){
                hidari_atom = &(Box->elements[i].atoms[j]);
                mass_scale += Box->elements[i].mass;
                Box->elements[i].atoms[j].type = hidari;
            }
            else if (current == middle_atom){
                naka_atom = &(Box->elements[i].atoms[j]);
                mass_scale += Box->elements[i].mass;
                Box->elements[i].atoms[j].type = naka;
                naka_mass = Box->elements[i].mass;
            }
            else if (current == after_atom){
                migi_atom = &(Box->elements[i].atoms[j]);
                mass_scale += Box->elements[i].mass;
                Box->elements[i].atoms[j].type = migi;
            }
            ++current;
        }
    }
    mass_scale /= naka_mass;
}

double *HMC_NVT_RC::adjust_angle(atom *move_atom, int type) const{
    Vector3<double> distance = move_atom->r - naka_atom->r;
    auto *sin_theta = new double [2];
    double curr_theta, new_theta, curr_phi, new_phi;
    switch (type)
    {
    case theta_1:
    case theta_2:
        sin_theta[0] = calculate_sin_theta(distance);
        curr_theta = acos(distance.z/distance.norm());
        new_theta = curr_theta + (Box->Random()-0.5) * theta_step;
        sin_theta[1] = sin(new_theta);
        move_atom->r.x = naka_atom->r.x + distance.x / sin_theta[0] * sin_theta[1];
        move_atom->r.y = naka_atom->r.y + distance.y / sin_theta[0] * sin_theta[1];
        move_atom->r.z = naka_atom->r.z + distance.norm() * cos(new_theta);
        return sin_theta;
        break;
    case phi_1:
    case phi_2:
        curr_theta = acos(distance.z/distance.norm());
        curr_phi = acos(distance.x/distance.norm()/sin(curr_theta));
        new_phi = curr_phi + (Box->Random()-0.5) * phi_step;
        move_atom->r.x = naka_atom->r.x + distance.x * cos(new_phi) / cos(curr_phi);
        move_atom->r.y = naka_atom->r.y + distance.y * sin(new_phi) / sin(curr_phi);
    default:
        return sin_theta;
        break;
    }
}

double HMC_NVT_RC::calculate_sin_theta(Vector3<double> vec){
    return sqrt(1-pow(vec.z/vec.norm(), 2));
}

HMC_NVT_DIS::HMC_NVT_DIS(box *Box, double (*calculate)(box *)){
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    set_constraint_atoms();
    if (!Consts.set_rc){
        Consts.reaction_coordinate = ::react_coord(correlation_atom, constraint_atom);
    }
    react_coord = Consts.reaction_coordinate;
    hardBoundary = new wall {Box};
}

void HMC_NVT_DIS::evolve(int &current, bool MC_step) const{
    if (MC_step)
        MC_evolve(current);
    else
        HMC_evolve(current);
}

void HMC_NVT_DIS::MC_evolve(int &current) const {
    double temp_energy = Box->potential_energy;
    Vector3<double> temp_pos = constraint_atom->r;
    Vector3<double> current_distance = constraint_atom->r - correlation_atom->r;
    double current_sin_theta = calculate_sin_theta(current_distance);
    double new_sin_theta = adjust_angle();
    if (hardBoundary->judge_accept()){
        Box->potential_energy = calculate(Box);
        if (Box->Random() > min(1.0, new_sin_theta/current_sin_theta*exp(-(Box->potential_energy-temp_energy)/(Box->k_B*Box->T)))){
            constraint_atom->r = temp_pos;
            Box->potential_energy = calculate(Box);
            current++;
        }
    }
    else {
        constraint_atom->r = temp_pos;
        current++;
    }
    Box->energy = Box->kinetic_energy + Box->potential_energy;
}

void HMC_NVT_DIS::HMC_evolve(int &current) const {
    double temp_energy = Box->energy;
    auto *temp_atoms = Box->new_atoms(Box->elements);
    correlation_atom->v *= (1/pow(mass_scale, 0.5)); // velocity is sampled for single atom, so it's necessary to scale for the rigid part.
    for (int step = 0 ; step < n_step; ++step) {
        for (int j = 0 ; j < Box->ntype ; ++j){
            element *curr_ele = &(Box->elements[j]);
            for (int k = 0 ; k < curr_ele->number ; ++k){
                if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                    if (curr_ele->atoms[k].type == hidari)
                        curr_ele->atoms[k].v += (correlation_atom->f+constraint_atom->f) * (delta_t*0.5) / (curr_ele->mass*mass_scale);
                    else if (curr_ele->atoms[k].type == normal)
                        curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t*0.5) / curr_ele->mass;
                }
            }
        }
        for (int j = 0 ; j < Box->ntype ; ++j){
            element *curr_ele = &(Box->elements[j]);
            for (int k = 0 ; k < curr_ele->number ; ++k){
                if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                    if (curr_ele->atoms[k].type != migi)
                        curr_ele->atoms[k].r += curr_ele->atoms[k].v * delta_t;
                    else
                        curr_ele->atoms[k].r += correlation_atom->v * delta_t;
                }
            }
        }
        Box->potential_energy = calculate(Box);
        for (int j = 0 ; j < Box->ntype ; ++j){
            element *curr_ele = &(Box->elements[j]);
            for (int k = 0 ; k < curr_ele->number ; ++k)
                if (curr_ele->atoms[k].move.x * curr_ele->atoms[k].move.y * curr_ele->atoms[k].move.z){
                    if (curr_ele->atoms[k].type == hidari)
                        curr_ele->atoms[k].v += (correlation_atom->f+constraint_atom->f) * (delta_t*0.5) / (curr_ele->mass*mass_scale);
                    else if (curr_ele->atoms[k].type == normal)
                        curr_ele->atoms[k].v += curr_ele->atoms[k].f * (delta_t*0.5) / curr_ele->mass;
                }
        }
    }
    correlation_atom->v *= pow(mass_scale, 0.5); // scale the velocity back to calculation kinetic energy
    Box->update_kinetic_energy();
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    if (Box->Random() > min(1.0, exp(-(Box->energy-temp_energy)/(Box->k_B*Box->T))) || !(hardBoundary->judge_accept())){
        Box->reload_r(temp_atoms);
        Box->potential_energy = calculate(Box);
        current++;
    }
    Box->init_velocities();
    Box->update_kinetic_energy();
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    delete [] temp_atoms;
}

void HMC_NVT_DIS::set_constraint_atoms(){
    int current = 0;
    int max_index = max(correlation_index, constraint_index);
    mass_scale = 0;
    double naka_mass;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            if (current > max_index)
                break;
            if (current == correlation_index){
                correlation_atom = &(Box->elements[i].atoms[j]);
                mass_scale += Box->elements[i].mass;
                Box->elements[i].atoms[j].type = hidari;
                naka_mass = Box->elements[i].mass;
            }
            else if (current == constraint_index){
                constraint_atom = &(Box->elements[i].atoms[j]);
                mass_scale += Box->elements[i].mass;
                Box->elements[i].atoms[j].type = migi;
            }
            ++current;
        }
    }
    mass_scale /= naka_mass;
}

double HMC_NVT_DIS::adjust_angle() const{
    Vector3<double> ball_move = sphere.ran_vector();
    constraint_atom->r += ball_move;
    Vector3<double> current_distance = constraint_atom->r - correlation_atom->r;
    double current_length = current_distance.norm();
    constraint_atom->r = correlation_atom->r + current_distance * react_coord / current_length;
    double sin_theta = calculate_sin_theta(current_distance);
    return sin_theta;
}