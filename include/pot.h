#ifndef CPIHMC_POT_H
#define CPIHMC_POT_H

#include "const.h"
#include "box.h"

namespace cpihmc
{
    class pot
    // parent class for potential energy and force calculator
    {
        public:
            explicit pot() = default;
            virtual ~pot() = default;
            virtual prec_t infer(box *) = 0;
    };
}

#endif