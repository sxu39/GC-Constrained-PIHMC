#ifndef CPIHMC_MD_NVE_H
#define CPIHMC_MD_NVE_H

#include "md.h"

namespace cpihmc
{
    class md_nve : public md
    {
        public:
            explicit md_nve(box * const Box, pot * const Pot, prec_t TimeStep):md(Box, Pot, TimeStep){}
            void evolve() const final;
    };
}

#endif