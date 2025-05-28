#ifndef CPIHMC_PIMC_H
#define CPIHMC_PIMC_H

#include "mc.h"

namespace cpihmc
{
    class pimc : public mc
    {
        private:
            const size_t ChangeBeadNum;
            std::vector<std::vector<vec3<prec_t>>> TempBeadCoords;
        private:
            void u2x(bead * const, const index_t) const;
        public:
            explicit pimc(const mc_mbr_pck &McMbrPck, const size_t ChangeBeadNum):mc(McMbrPck), ChangeBeadNum(ChangeBeadNum), TempBeadCoords(Box->Beads.size()){}
            ~pimc() = default;
            const bool_t mc_evolve() final;
    };
}

#endif