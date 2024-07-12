//
// Created by Jin Bin on 2022/02/10.
//

#ifndef MD_MC_FIRST_PRINCIPLE_H
#define MD_MC_FIRST_PRINCIPLE_H

#include "box.h"
#include "cell.h"
#include <cstdlib>
#include "potential.h"

extern string format;

class vasp : public potential{
public:
    static double calculate(box *Box);
};

class abacus : public potential{
public:
    static double calculate(box *Box);
};

#endif //MD_MC_FIRST_PRINCIPLE_H
