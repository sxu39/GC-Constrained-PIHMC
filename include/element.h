//
// Created by Jin Bin on 2022/02/02.
//

#ifndef MD_MC_ELEMENT_H
#define MD_MC_ELEMENT_H

#include "atom.h"

class element{
public:
    string label;
    int number;
    double mass;
    atom *atoms;
    element();
    ~element();
};

#endif