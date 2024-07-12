//
// Created by Jin Bin on 2022/02/10.
//

#ifndef CPIHMC_DP_POTENTIAL_H
#define CPIHMC_DP_POTENTIAL_H

#include "box.h"
#include "cell.h"
#include "deepmd/DeepPot.h"
#include <cstdlib>
#include "potential.h"

typedef double ValueType;

extern deepmd::DeepPot *dp;
extern deepmd::DeepPot canonical_dp;
extern deepmd::DeepPot grand_canonical_dp;

class deep_potential : public potential {
public:
    static double calculate(box *Box);
    friend class model_devi;
private:
    static void provide_structures(vector<ValueType> &, vector<int> &, vector<ValueType> &, box *);
    static void single_point_calculate(vector<double> &, vector<ValueType> &, const vector<ValueType> &, const vector<int> &, const vector<ValueType> &);
    static void single_point_calculate(vector<double> &, vector<ValueType> &, const vector<ValueType> &, const vector<int> &, const vector<ValueType> &, const vector<ValueType> &);
};


#endif //CPIHMC_DP_POTENTIAL_H
