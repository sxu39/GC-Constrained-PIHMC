#ifndef CPIHMC_HMC_H
#define CPIHMC_HMC_H

#include "rand.h"

namespace cpihmc
{
    class hmc
    {
        protected:
            const prec_t Temp;
            const prec_t TimeStep;
            const prec_t MassScal;
            rand &Rand;
        public:
            explicit hmc(const prec_t Temp, const prec_t TimeStep, const prec_t MassScal, rand &Rand):Temp(Temp), TimeStep(TimeStep), MassScal(MassScal), Rand(Rand){}
            virtual ~hmc() = default;
            virtual void random_velocities() = 0;
            virtual void hmc_evolve_forth_before() = 0;
            virtual void hmc_evolve_forth_after() = 0;
            virtual void hmc_evolve_back() = 0;
    };
}

#endif