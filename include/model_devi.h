//
// Created by Jin Bin on 2023/10/21.
//

#ifndef CPIHMC_MODEL_DEVI_H
#define CPIHMC_MODEL_DEVI_H

#include "box.h"
#include "ml_potential.h"

class model_devi {
public:
    model_devi() = default;
    model_devi(box *);
    void init(box *);
    void compute_model_devi(int);
    ~model_devi();
private:
    box *Box;
    int num_models;
    vector<deepmd::DeepPot> model_devi_models;
    ofstream model_devi_output;
    static void ana_st(double &, double &,double &, const vector<double> &, const int &);
    void compute_force_virials(vector<vector<ValueType>> &, vector<vector<ValueType>> &, 
    const vector<ValueType> &, const vector<int> &, const vector<ValueType> &);
    void compute_force_virials(vector<vector<ValueType>> &, vector<vector<ValueType>> &, 
    const vector<ValueType> &, const vector<int> &, const vector<ValueType> &, const vector<ValueType> &);
    void cope_with_force(double &, double &, double &, const vector<vector<double>> &);
    void cope_with_virial(double &, double &, double &, const vector<vector<double>> &);
    template <typename VALUETYPE>
    void compute_avg(vector<VALUETYPE> &, const vector<vector<VALUETYPE>> &);
    template <typename VALUETYPE>
    void compute_std(vector<VALUETYPE> &, const vector<VALUETYPE> &, const vector<vector<VALUETYPE>> &, 
    const int &);
};

#endif // CPIHMC_MODEL_DEVI_H