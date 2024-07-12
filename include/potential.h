//
// Created by Jin Bin on 2022/09/09.
//

#ifndef CPIHMC_POTENTIAL_H
#define CPIHMC_POTENTIAL_H

class potential
// parent class for potential energy and force calculator
{
public:
    static double eV2energy; // energy transformation coefficient, from eV
    static double Ry_d_Bohr2force; // force transformation coefficient, from Ry/Bohr
    static double eV_d_A2force; // force transformation coefficient, from eV/Ang
    static double Bohr2A; // length transformation coefficient, from Bohr to Ang
};

#endif //CPIHMC_POTENTIAL_H