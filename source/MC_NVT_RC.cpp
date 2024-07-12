//
// Created by Jin Bin on 2021/11/03.
//

#include "MC_NVT_RC.h"
#include "gen_func.h"

MC_NVT_RC::MC_NVT_RC(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    number_list = new int [Box->ntype+1];
    number_list[0] = 0;
    int current = 0;
    for (int i = 0 ; i < Box->ntype ; ++i){
        current += Box->elements[i].number;
        number_list[i+1] = current;
    }
    set_constraint_atoms();
    if (!Consts.set_rc){
        Consts.reaction_coordinate = ::react_coord(hidari_atom, naka_atom, migi_atom);
    }
    reaction_coordinate = Consts.reaction_coordinate;
}

MC_NVT_RC::~MC_NVT_RC(){
    delete [] number_list;
}

void MC_NVT_RC::evolve(int &current) const {
    double current_energy = Box->potential_energy;
    atom *change_atom = change_atom_pointer();
    Vector3<double> current_pos = change_atom->r;
    double factor = 1;
    if (change_atom->type == naka)
        scale(factor);
    else{
        Vector3<double> trial_move_direction = trial_move.ran_vector();
        change_atom->r += trial_move_direction;
    }
    if (change_atom->type % 2){
        double curr_dis = (current_pos - naka_atom->r).norm();
        scale(naka_atom, change_atom, curr_dis);
    }
    Box->potential_energy = calculate(Box);
    if (Box->Random() > min(1.0, factor*exp(-(Box->potential_energy-current_energy)/(Box->k_B*Box->T)))){
        change_atom->r = current_pos;
        Box->potential_energy = calculate(Box);
        current++;
    }
    // delete [] current_atoms;
}

template<class T>
T MC_NVT_RC::min(T x, T y) {
    return x<y?x:y;
}

atom *MC_NVT_RC::change_atom_pointer() const {
    atom *change_atom;
    int atom_index;
    do{
        atom_index = floor(Box->Random()*Box->N_atoms);
        for (int i = 0 ; i < Box->ntype ; ++i){
            if (atom_index < number_list[i+1]){
                change_atom = &(Box->elements[i].atoms[atom_index-number_list[i]]);
                break;
            }
        }
    } while (!(change_atom->move.x * change_atom->move.y * change_atom->move.z));
    return change_atom;
}

void MC_NVT_RC::set_constraint_atoms(){
    int current = 0;
    int max_index = max(left_atom, max(middle_atom, after_atom));
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            if (current > max_index)
                break;
            if (current == left_atom){
                hidari_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = hidari;
            }
            else if (current == middle_atom){
                naka_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = naka;
            }
            else if (current == after_atom){
                migi_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = migi;
            }
            ++current;
        }
    }
}

void MC_NVT_RC::scale(double &factor) const{
    Vector3<double> axis = migi_atom->r - hidari_atom->r; // (a, b, c)
    Matrix3 rotation;
    double norm3 = axis.norm(); // \sqrt{a^2+b^2+c^2}
    if (axis.y){
        double norm2 = pow(pow(axis.y, 2)+pow(axis.z, 2), 0.5); // \sqrt{b^2+c^2}
        rotation.e11 = norm2 / norm3;
        rotation.e12 = -axis.x * axis.y / (norm2 * norm3);
        rotation.e13 = -axis.x * axis.z / (norm2 * norm3);
        rotation.e22 = axis.z / norm2;
        rotation.e23 = -axis.y / norm2;
        rotation.e31 = axis.x / norm3;
        rotation.e32 = axis.y / norm3;
        rotation.e33 = axis.z / norm3;
    }
    else{
        double norm2 = pow(pow(axis.x, 2)+pow(axis.z, 2), 0.5); // \sqrt{a^2+c^2}
        rotation.e11 = axis.z / norm2;
        rotation.e13 = -axis.x / norm2;
        rotation.e31 = axis.x / norm2;
        rotation.e33 = axis.z / norm2;
    }
    Vector3<double> temp_hidari = rotation * hidari_atom->r; // coordinate after rotation
    Vector3<double> temp_naka = rotation * naka_atom->r;
    Vector3<double> temp_migi = rotation * migi_atom->r;
    temp_naka -= temp_migi; // the vector from migi atom to naka atom, which is used to calculation theta and phi
    double c = norm3 / 2; // constant for conic section
    double a = reaction_coordinate / 2;
    double curr_rho = (naka_atom->r-migi_atom->r).norm();
    double curr_theta = acos(temp_naka.z / curr_rho);
    double new_theta = curr_theta + (Box->Random()-0.5) * theta_step;
    double new_rho = (pow(c, 2) - pow(a, 2))/(a-c*cos(new_theta));
    factor = (new_rho+reaction_coordinate) / (curr_rho+reaction_coordinate) * pow(new_rho/curr_rho, 3) * sin(new_theta) / sin(curr_theta);
    double delta_phi = (Box->Random()-0.5) * phi_step;
    double naka_norm2 = pow(pow(temp_naka.x, 2)+pow(temp_naka.y, 2), 0.5);
    double cos_phi = temp_naka.x / naka_norm2;
    double sin_phi = temp_naka.y / naka_norm2;
    Vector3<double> new_naka;
    new_naka.x = temp_migi.x + new_rho * sin(new_theta) * (cos_phi * cos(delta_phi) - sin_phi * sin(delta_phi));
    new_naka.y = temp_migi.y + new_rho * sin(new_theta) * (sin_phi * cos(delta_phi) + cos_phi * sin(delta_phi));
    new_naka.z = temp_migi.z + new_rho * cos(new_theta);
    naka_atom->r = (rotation.Inverse()) * new_naka;
}

void MC_NVT_RC::scale(atom *correlation_atom, atom *constraint_atom, double length) const {
    Vector3<double> current_distance = constraint_atom->r - correlation_atom->r;
    double current_length = current_distance.norm();
    constraint_atom->r = correlation_atom->r + current_distance * length / current_length;
}

double MC_NVT_RC::free_energy_estimator_left() const{
    double free_energy = 0;
    Vector3<double> direction = hidari_atom->r - naka_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy -= direction * hidari_atom->f;
    free_energy -= 2*Box->k_B*Box->T * length_inv;
    return free_energy;
}

double MC_NVT_RC::free_energy_estimator_right() const{
    double free_energy = 0;
    Vector3<double> direction = migi_atom->r - naka_atom->r;
    double length_inv = 1/direction.norm();
    direction *= length_inv;
    free_energy -= direction * migi_atom->f;
    free_energy -= 2*Box->k_B*Box->T * length_inv;
    return free_energy;
}

MC_NVT_DIS::MC_NVT_DIS(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    Box->potential_energy = calculate(Box);
    number_list = new int [Box->ntype+1];
    number_list[0] = 0;
    int current = 0;
    for (int i = 0 ; i < Box->ntype ; ++i){
        current += Box->elements[i].number;
        number_list[i+1] = current;
    }
    set_constraint_atoms();
    if (!Consts.set_rc){
        Consts.reaction_coordinate = ::react_coord(correlation_atom, constraint_atom);
    }
    reaction_coordinate = Consts.reaction_coordinate;
}

void MC_NVT_DIS::evolve(int &current) const{
    double current_energy = Box->potential_energy;
    atom *change_atom = change_atom_pointer();
    Vector3<double> current_pos = change_atom->r;
    double factor = 1;
    Vector3<double> trial_move_direction = trial_move.ran_vector();
    change_atom->r += trial_move_direction;
    if (change_atom->type == hidari){
        double curr_dis = (current_pos - constraint_atom->r).norm();
        scale(constraint_atom, change_atom, curr_dis);
    }
    else if (change_atom->type == migi){
        double curr_dis = (current_pos - correlation_atom->r).norm();
        scale(correlation_atom, change_atom, curr_dis);
    }
    Box->potential_energy = calculate(Box);
    if (Box->Random() > min(1.0, factor*exp(-(Box->potential_energy-current_energy)/(Box->k_B*Box->T)))){
        change_atom->r = current_pos;
        Box->potential_energy = calculate(Box);
        current++;
    }
    // delete [] current_atoms;
}

void MC_NVT_DIS::set_constraint_atoms(){
    int current = 0;
    int max_index = max(correlation_index, constraint_index);
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            if (current > max_index)
                break;
            if (current == correlation_index){
                correlation_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = hidari;
            }
            else if (current == constraint_index){
                constraint_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = migi;
            }
            ++current;
        }
    }
}