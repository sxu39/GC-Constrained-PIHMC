#include "model_devi.h"

model_devi::model_devi(box *Box):Box(Box), num_models(Consts.model_devi_DP_models.size()), model_devi_output(Consts.model_devi_file){
    model_devi_models.resize(num_models);
    for (int i = 0 ; i < num_models ; ++i)
        model_devi_models[i].init(Consts.model_devi_DP_models[i]);
    model_devi_output << scientific;
    model_devi_output << "#" << setw(12 - 1) << "step" << setw(18 + 1) << "max_devi_v"
        << setw(18 + 1) << "min_devi_v" << setw(18 + 1) << "avg_devi_v"
        << setw(18 + 1) << "max_devi_f" << setw(18 + 1) << "min_devi_f"
        << setw(18 + 1) << "avg_devi_f" << endl;
}

void model_devi::init(box *Box){
    this->Box = Box;
    num_models = Consts.model_devi_DP_models.size();
    model_devi_models.resize(num_models);
    model_devi_output.open(Consts.model_devi_file);
    for (int i = 0 ; i < num_models ; ++i)
        model_devi_models[i].init(Consts.model_devi_DP_models[i]);
    model_devi_output << scientific;
    model_devi_output << "#" << setw(12 - 1) << "step" << setw(18 + 1) << "max_devi_v"
        << setw(18 + 1) << "min_devi_v" << setw(18 + 1) << "avg_devi_v"
        << setw(18 + 1) << "max_devi_f" << setw(18 + 1) << "min_devi_f"
        << setw(18 + 1) << "avg_devi_f" << endl;
}

void model_devi::compute_model_devi(int n_step) {
    // collect all beads structures in the format required for DP inference
    vector<ValueType> coords(Consts.P*Box->N_atoms*3);
    vector<ValueType> cells;
    vector<int> atype(Box->N_atoms);
    deep_potential::provide_structures(coords, atype, cells, Box);

    vector<vector<ValueType>> all_force;
    vector<vector<ValueType>> all_virial;

    if (strcmp(Consts.GC_type.c_str(), "Continuous") == 0){
        vector<ValueType> electron_numbers(Consts.P, *Box->electron_number);
        compute_force_virials(all_force, all_virial, coords, atype, cells, electron_numbers);
    }
    else if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0)
        compute_force_virials(all_force, all_virial, coords, atype, cells);

    double all_f_min = numeric_limits<double>::max(), all_f_max = 0, all_f_avg = 0;
    cope_with_force(all_f_max, all_f_min, all_f_avg, all_force);

    for (int kk = 0; kk < num_models; ++kk) {
        for (int ii = 0; ii < 9 * Consts.P; ++ii)
            all_virial[kk][ii] /= Box->N_atoms;
    }
    double all_v_min = numeric_limits<double>::max(), all_v_max = 0, all_v_avg = 0;
    cope_with_virial(all_v_max, all_v_min, all_v_avg, all_virial);

    model_devi_output << setw(12) << n_step << " " << setw(18) << all_v_max
        << " " << setw(18) << all_v_min << " " << setw(18) << all_v_avg
        << " " << setw(18) << all_f_max << " " << setw(18) << all_f_min
        << " " << setw(18) << all_f_avg << endl;
}

model_devi::~model_devi(){
    if (model_devi_output.is_open())
        model_devi_output.close();
}

void model_devi::ana_st(double &max, double &min, double &avg, const vector<double> &vec, const int &nloc) {
    if (nloc == 0) return;
    max = vec[0];
    min = vec[0];
    avg = vec[0];
    for (unsigned ii = 1; ii < nloc; ++ii) {
        if (vec[ii] > max) {
        max = vec[ii];
        }
        if (vec[ii] < min) {
        min = vec[ii];
        }
        avg += vec[ii];
    }
    avg /= nloc;
}

void model_devi::compute_force_virials(vector<vector<ValueType>> &forces, vector<vector<ValueType>> &virials, 
const vector<ValueType> &coords, const vector<int> &atype, const vector<ValueType> &cells){
    vector<double> energies;
    forces.resize(num_models);
    virials.resize(num_models);
    for (int i = 0 ; i < num_models ; ++i)
        model_devi_models[i].compute(energies, forces[i], virials[i], coords, atype, cells);
}

void model_devi::compute_force_virials(vector<vector<ValueType>> &forces, vector<vector<ValueType>> &virials, 
const vector<ValueType> &coords, const vector<int> &atype, const vector<ValueType> &cells, const vector<ValueType> &electron_numbers){
    vector<double> energies;
    forces.resize(num_models);
    virials.resize(num_models);
    for (int i = 0 ; i < num_models ; ++i)
        model_devi_models[i].compute(energies, forces[i], virials[i], coords, atype, cells, electron_numbers);
}

void model_devi::cope_with_force(double &all_f_max, double &all_f_min, double &all_f_avg, const vector<vector<double>> &all_force){
    vector<double> tmp_avg_f;
    vector<double> std_f;
    compute_avg(tmp_avg_f, all_force);
    compute_std(std_f, tmp_avg_f, all_force, 3);
    ana_st(all_f_max, all_f_min, all_f_avg, std_f, Box->N_atoms*Consts.P);
    all_f_max *= potential::eV_d_A2force;
    all_f_min *= potential::eV_d_A2force;
    all_f_avg *= potential::eV_d_A2force;
}

void model_devi::cope_with_virial(double &all_v_max, double &all_v_min, double &all_v_avg, const vector<vector<double>> &all_virial){
    vector<double> avg_virial;
    vector<double> std_virial;
    compute_avg(avg_virial, all_virial);
    compute_std(std_virial, avg_virial, all_virial, 1);
    for (int ii = 0; ii < 9 * Consts.P; ++ii) {
        if (std_virial[ii] > all_v_max) {
            all_v_max = std_virial[ii];
        }
        if (std_virial[ii] < all_v_min) {
            all_v_min = std_virial[ii];
        }
        all_v_avg += std_virial[ii] * std_virial[ii];
    }
    all_v_avg = sqrt(all_v_avg / (9*Consts.P));
    all_v_max *= potential::eV2energy;
    all_v_min *= potential::eV2energy;
    all_v_avg *= potential::eV2energy;
}

template <typename VALUETYPE>
void model_devi::compute_avg(vector<VALUETYPE> &avg, const vector<vector<VALUETYPE>> &xx) {
    assert(xx.size() == num_models);
    if (num_models == 0) {
        return;
    }

    avg.resize(xx[0].size());
    fill(avg.begin(), avg.end(), VALUETYPE(0.));

    for (unsigned ii = 0; ii < num_models; ++ii) {
        for (unsigned jj = 0; jj < avg.size(); ++jj) {
            avg[jj] += xx[ii][jj];
        }
    }

    for (unsigned jj = 0; jj < avg.size(); ++jj) {
        avg[jj] /= VALUETYPE(num_models);
    }
}

template void model_devi::compute_avg<double>(vector<double> &avg, const vector<vector<double>> &xx);

template <typename VALUETYPE>
void model_devi::compute_std(vector<VALUETYPE>& std, const vector<VALUETYPE>& avg, 
const vector<vector<VALUETYPE>>& xx, const int &stride) {
    assert(xx.size() == num_models);
    if (num_models == 0){
        return;
    }

    unsigned ndof = avg.size();
    unsigned nloc = ndof / stride;
    assert(nloc * stride == ndof);

    std.resize(nloc);
    fill(std.begin(), std.end(), VALUETYPE(0.));

    for (unsigned ii = 0; ii < num_models; ++ii) {
        for (unsigned jj = 0; jj < nloc; ++jj) {
            const VALUETYPE* tmp_f = &(xx[ii][jj * stride]);
            const VALUETYPE* tmp_avg = &(avg[jj * stride]);
            for (unsigned dd = 0; dd < stride; ++dd) {
                VALUETYPE vdiff = tmp_f[dd] - tmp_avg[dd];
                std[jj] += vdiff * vdiff;
            }
        }
    }

    for (unsigned jj = 0; jj < nloc; ++jj) {
        std[jj] = sqrt(std[jj] / VALUETYPE(num_models));
    }
}

template void model_devi::compute_std<double>(vector<double>& std, const vector<double>& avg,
const vector<vector<double>>& xx, const int &stride);