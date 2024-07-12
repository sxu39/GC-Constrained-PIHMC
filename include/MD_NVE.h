//
// Created by Jin Bin on 2022/02/21.
//

#ifndef CPIHMC_MD_NVE_H
#define CPIHMC_MD_NVE_H

#include "box.h"

class MD_NVE {
public:
    box *Box;
    double (*calculate)(box *);
    explicit MD_NVE(box *Box, double (*calculate)(box *));
    void evolve(double delta_t) const;
};


#endif //CPIHMC_MD_NVE_H
