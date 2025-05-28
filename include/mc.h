#ifndef CPIHMC_MC_H
#define CPIHMC_MC_H

#include "pot.h"
#include "temp.h"
#include "hard_bdry.h"

namespace cpihmc
{
    struct mc_mbr_pck
    {
        box * const Box;
        pot * const Pot;
        temp * const Temp;
        hard_bdry * const HardBdry;
        rand &Rand;
        explicit mc_mbr_pck(box * const Box, pot * const Pot, temp * const Temp, hard_bdry * const HardBdry, rand &Rand):Box(Box), Pot(Pot), Temp(Temp), HardBdry(HardBdry), Rand(Rand){}
    };

    class mc
    {
        protected:
            box * const Box;
            pot * const Pot;
            temp * const Temp;
            hard_bdry * const HardBdry;
            rand &Rand;
        public:
            explicit mc(const mc_mbr_pck &McMbrPck):Box(McMbrPck.Box), Pot(McMbrPck.Pot), Temp(McMbrPck.Temp),HardBdry(McMbrPck.HardBdry), Rand(McMbrPck.Rand){}
            virtual ~mc() = default;
            virtual const bool_t mc_evolve() = 0;
    };
}

#endif