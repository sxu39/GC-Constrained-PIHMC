#ifndef CPIHMC_HARD_BDRY_H
#define CPIHMC_HARD_BDRY_H

#include "box.h"

namespace cpihmc
{
    class hard_bdry
    {
        protected:
            box * const Box; // simulated system
        public:
            explicit hard_bdry(box *Box):Box(Box){}
            virtual ~hard_bdry() = default;
            virtual const bool_t judge_accept() const = 0;
    };
}

#endif