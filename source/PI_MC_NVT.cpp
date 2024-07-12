//
// Created by Jin Bin on 2022/02/21.
//

#include "PI_MC_NVT.h"
#include "gen_func.h"

#define hbar 0.0226974

PI_MC_NVT::PI_MC_NVT(box *Box, double (*calculate)(box *)) {
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    // calculate system energy and force for the first time
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
    react_coord = Consts.reaction_coordinate;
}

PI_MC_NVT::~PI_MC_NVT(){
    // string command = "rm " + Consts.work_path + "_* -r";
    // system(command.c_str());
    delete [] number_list;
}

void PI_MC_NVT::evolve(int &current, bool centroid) const{
    if (centroid)
        centroid_evolve(current);
    else
        internal_evolve(current);
}

void PI_MC_NVT::internal_evolve(int &current) const{
    double temp_energy = Box->potential_energy;
    double temp_quantum_kinetic_energy = Box->quantum_kinetic_energy;
    int bead_num = Box->beads.size();
    vector<vector<Vector3<double> > > curr_bead_pos(bead_num);
    for (int i = 0 ; i < bead_num ; ++i)
        curr_bead_pos[i] = Box->beads[i].r;

    vector<Vector3<double> > curr_force(Box->N_atoms); // to save forces of all atom in current structure
    int curr_index = 0;
    for (int i = 0 ; i < Box->ntype ; ++i){
        for (int j = 0 ; j < Box->elements[i].number ; ++j){
            curr_force[curr_index] = Box->elements[i].atoms[j].f;
            ++curr_index;
        }
    }
    
    int change_bead = floor(Box->Random()*P);
    for (int i = 0 ; i < bead_num ; ++i){
        u2x(&(Box->beads[i]), change_bead);
        back_centroid(&(Box->beads[i]));
    }
    Box->potential_energy = calculate(Box);
    if (Box->Random() > min(1.0, exp(-(Box->potential_energy-temp_energy)/Box->k_B/Box->T))){
        for (int i = 0 ; i < bead_num ; ++i)
            Box->beads[i].r = curr_bead_pos[i];
        Box->potential_energy = temp_energy;
        Box->quantum_kinetic_energy = temp_quantum_kinetic_energy;
        curr_index = 0;
        for (int i = 0 ; i < Box->ntype ; ++i){
            for (int j = 0 ; j < Box->elements[i].number ; ++j){
                Box->elements[i].atoms[j].f = curr_force[curr_index];
                ++curr_index;
            }
        }
        current++;
    }
    Box->energy = Box->kinetic_energy + Box->potential_energy;
}

void PI_MC_NVT::centroid_evolve(int &current) const{
    double current_energy = Box->potential_energy;
    double temp_quantum_kinetic_energy = Box->quantum_kinetic_energy;
    atom *change_atom = change_atom_pointer();
    Vector3<double> current_pos = change_atom->r;
    double factor = 1;
    if (change_atom->type == naka){
        scale(factor);
        back_centroid(change_atom->pbead);
    }
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

        if (change_atom->type == naka)
            back_centroid(change_atom->pbead); // only for one bead atom situation
        
        Box->potential_energy = calculate(Box);
        Box->quantum_kinetic_energy = temp_quantum_kinetic_energy;
        current++;
    }
}

atom *PI_MC_NVT::change_atom_pointer() const {
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

void PI_MC_NVT::u2x(bead *bead, int change_bead) const{
    double PI_mass = bead->mass;
    for (int i = J ; i > 0 ; i--){
        Vector3<double> u;
        u.x = Box->RanGaussian(PI_mass*(i+1)*P*pow(Box->k_B*Box->T, 2)/(i*pow(hbar, 2)));
        u.y = Box->RanGaussian(PI_mass*(i+1)*P*pow(Box->k_B*Box->T, 2)/(i*pow(hbar, 2)));
        u.z = Box->RanGaussian(PI_mass*(i+1)*P*pow(Box->k_B*Box->T, 2)/(i*pow(hbar, 2)));
        // temp[i] = u + (float(i) * temp[i + 1] + temp[0]) / (float(i) + 1); // transform
        bead->r[(change_bead+i)%P] = u + ((double)i*bead->r[(change_bead+i+1)%P]+bead->r[change_bead%P]) * (1/double(i+1));
    }
}

void PI_MC_NVT::back_centroid(bead *bead) const {
    Vector3<double> displacement = bead->bead_atom->r;
    double P_inv = 1/(double)P;
    for (int i = 0 ; i < P ; ++i)
        displacement -= bead->r[i] * P_inv;
    for (int i = 0 ; i < P ; ++i)
        bead->r[i] += displacement;
}

void PI_MC_NVT::set_constraint_atoms(){
    int current = 0;
    int max_index = max(left_atom, max(middle_atom, after_atom));
    double naka_mass;
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
                naka_mass = Box->elements[i].mass;
            }
            else if (current == after_atom){
                migi_atom = &(Box->elements[i].atoms[j]);
                Box->elements[i].atoms[j].type = migi;
            }
            ++current;
        }
    }
}

void PI_MC_NVT::scale(double &factor) const{
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
    double a = react_coord/ 2;
    double curr_rho = (naka_atom->r-migi_atom->r).norm();
    double curr_theta = acos(temp_naka.z / curr_rho);
    double new_theta = curr_theta + (Box->Random()-0.5) * theta_step;
    double new_rho = (pow(c, 2) - pow(a, 2))/(a-c*cos(new_theta));
    factor = (new_rho+react_coord) / (curr_rho+react_coord) * pow(new_rho/curr_rho, 3) * sin(new_theta) / sin(curr_theta);
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

void PI_MC_NVT::scale(atom *correlation_atom, atom *constraint_atom, double length) const {
    Vector3<double> current_distance = constraint_atom->r - correlation_atom->r;
    double current_length = current_distance.norm();
    constraint_atom->r = correlation_atom->r + current_distance * length / current_length;
}
