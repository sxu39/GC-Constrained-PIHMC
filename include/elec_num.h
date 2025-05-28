#ifndef CPIHMC_ELEC_NUM_H
#define CPIHMC_ELEC_NUM_H

#include "mc.h"

namespace cpihmc
{
    class elec_num : public mc
    {
        private:
            const prec_t Mu;
            const std::pair<prec_t, prec_t> ElecNumRange;
            const prec_t ElecNumWidth;
        public:
            explicit elec_num(const mc_mbr_pck &McMbrPck, const prec_t Mu, const std::pair<prec_t, prec_t> &ElecNumRange, const prec_t ElecNumWidth):Mu(Mu*eV2Energy), ElecNumRange(ElecNumRange), ElecNumWidth(ElecNumWidth), mc(McMbrPck){}
            ~elec_num() = default;
            const bool_t mc_evolve() final;
    };
}

#endif