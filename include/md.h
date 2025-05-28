#ifndef CPIHMC_MD_H
#define CPIHMC_MD_H

#include "pot.h"

namespace cpihmc
{
    class md
    {
        protected:
            box * const Box;
            pot * const Pot;
            prec_t TimeStep;
        public:
            explicit md(box * const Box, pot * const Pot, prec_t TimeStep):Box(Box), Pot(Pot), TimeStep(TimeStep){}
            virtual ~md() = default;
            virtual void evolve() const = 0;
    };
}

#endif