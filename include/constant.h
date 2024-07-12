//
// Created by Jin Bin on 2021/10/03.
//

#ifndef MD_MC_CONSTANT_H
#define MD_MC_CONSTANT_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cassert>
using namespace std;

template <class T>
void OUTP(ofstream &ofs, const string &name, const T &a, const string &explanation="")
{
    ofs << setw(20) << name << a << " #" << explanation << endl;
}

class constant{
public:
    // General
    string atom_file;
    string work_path;
    string out_file;
    string DP_model;
    string calculation_type;
    string potential_type;
    int steps;
    int dump_steps;
    vector<string> output_terms;
    unordered_map<string, int> quantity_digits;
    // System
    double T;
    int ntype;
    vector<string> element_type;
    vector<int> element_index;
    string RC_type;
    int *constraint_atom;
    double reaction_coordinate;
    vector<string> wall;
    vector<double> wall_pos;
    int N_MAX;
    // Grand Cononical Ensemble
    double Ne_proportion;
    string DP_model_prefix;
    pair<double, double> Ne_range;
    double mu;
    double delta_Ne;
    string GC_type;
    // Model Deviation
    vector<string> model_devi_DP_models;
    int model_devi_interval;
    string model_devi_file;
    // Hamiltonian Monte Carlo
    double hmc_proportion;
    int N_step;
    double length_step;
    double M_scaling;
    // Hamiltonian Monte Carlo / Molecular Dynamics
    double Delta_t;
    // Monte Carlo
    double radial;
    // (Hamiltonian) Monte Carlo
    double theta_step;
    double phi_step;
    // Path Integral
    string initial_beads;
    string beads_file;
    int P;
    int J;
    int bead_num;
    vector<int> bead_index;
    constant();
    constant(const string &input);
    ~constant();
    void Read(const string &input);
    void Print(const string &output) const;
    static void strtolower(const char *sa, char *sb);
    bool set_rc;
private:
    template <class T>
    static void read_value(ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
    }
    template <class T>
    static void read_value(ifstream &ifs, vector<T> &var_vec)
    {
        T temp;
        while (ifs.peek() != '\n'){
            if (ifs.peek() == ' ' || ifs.peek() == '\t')
                ifs.get();
            else if (ifs.peek() == '#'){
                ifs.ignore(150, '\n');
                break;
            }
            else {
                ifs >> temp;
                var_vec.push_back(temp);
            }
        }
        if (ifs.peek() == '\n')
            ifs.get(); // get rid of '\n'
    }
};

#endif //MD_MC_CONSTANT_H
