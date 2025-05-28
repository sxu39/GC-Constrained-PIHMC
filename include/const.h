#ifndef CPIHMC_CONST_H
#define CPIHMC_CONST_H

#include "type.h"

namespace cpihmc
{
    const prec_t k_B = 2.96915e-6;
    const prec_t hbar = 0.0226791;
    const prec_t eV2Energy = 0.0344556; // energy transformation coefficient, from eV
    const prec_t Ry_d_Bohr2Force = 0.468792; // force transformation coefficient, from Ry/Bohr
    const prec_t eV_d_Ang2Force = 0.0182331; // force transformation coefficient, from eV/Ang
    const prec_t Bohr2Ang = 0.529177; // length transformation coefficient, from Bohr to Ang
}

#endif