//
// Created by Jin Bin on 2021/11/08.
//

#include "box.h"

// constant Consts; //NOLINT

box::box(element *Elements, Matrix3 Lattice, int N_atoms, double *electron_number, long long Seed) : elements(Elements), lattice_vector(Lattice), N_atoms(N_atoms),
electron_number(electron_number), seed(Seed){
//    N_atoms = int(atoms.size());
    init_velocities();
    update_kinetic_energy();
    if (Consts.element_index.size()){
        for (int i = 0 ; i < ntype ; ++i)
            element_type.insert(pair<string, int>(Consts.element_type[i], Consts.element_index[i]));
    }
    else {
        for (int i = 0 ; i < ntype ; ++i)
            element_type.insert(pair<string, int>(Consts.element_type[i], i));
    }
    lattice_vector_inverse = lattice_vector.Inverse();
}

void box::init_velocities() {
    for (int i = 0 ; i < ntype ; ++i){
        double mass = elements[i].mass;
        for (int j = 0 ; j < elements[i].number ; ++j){
            if (elements[i].atoms[j].move.x)
                elements[i].atoms[j].v.x = RanGaussian(mass);
            if (elements[i].atoms[j].move.y)
                elements[i].atoms[j].v.y = RanGaussian(mass);
            if (elements[i].atoms[j].move.z)
                elements[i].atoms[j].v.z = RanGaussian(mass);
        }
    }
}

void box::update_kinetic_energy() {
    double sum = 0;
    for (int i = 0 ; i < ntype ; ++i){
        double mass = elements[i].mass;
        for (int j = 0 ; j < elements[i].number ; ++j){
            if (elements[i].atoms[j].type == normal || elements[i].atoms[j].type == naka){
                sum += elements[i].atoms[j].v * elements[i].atoms[j].v * 0.5 * mass;
            }
        }
    }
    kinetic_energy = sum;
}

double box::Random() const{
    static default_random_engine eng(seed);
    static uniform_real_distribution<double> random(0.0, 1.0);
    return random(eng);
}

double *box::new_atoms(element *Elements) const {
    auto *temp_atoms = new double [3*N_atoms];
    int current = 0;
    for (int i = 0 ; i < ntype ; ++i) {
        for (int j = 0 ; j < elements[i].number ; ++j){
            temp_atoms[3*(current+j)] = elements[i].atoms[j].r.x;
            temp_atoms[3*(current+j)+1] = elements[i].atoms[j].r.y;
            temp_atoms[3*(current+j)+2] = elements[i].atoms[j].r.z;
        }
        current += elements[i].number;
    }
    return temp_atoms;
}

void box::reload_r(const double *vec) const {
    int current = 0;
    for (int i = 0 ; i < ntype ; ++i){
        for (int j = 0 ; j < elements[i].number ; ++j){
            elements[i].atoms[j].r.x = vec[3*(current+j)];
            elements[i].atoms[j].r.y = vec[3*(current+j)+1];
            elements[i].atoms[j].r.z = vec[3*(current+j)+2];
        }
        current += elements[i].number;
    }
}

void box::regular() const {
    for (int i = 0 ; i < ntype ; ++i){
        for (int j = 0 ; j < elements[i].number ; ++j)
            elements[i].atoms[j].r = ((elements[i].atoms[j].r * lattice_vector_inverse) % 1.0) * lattice_vector;
    }
}

double box::RanGaussian(double mass) const {
    double random_Gaussian = sqrt(T*k_B/mass) * std::sqrt(-2*std::log(Random()))*cos(2*PI*Random());
    // if the random number of Gaussian distribution is inf, than calculate one more time
    while (std::isinf(random_Gaussian)){
        random_Gaussian = sqrt(T*k_B/mass) * std::sqrt(-2*std::log(Random()))*cos(2*PI*Random());
    }
    return random_Gaussian;
}
