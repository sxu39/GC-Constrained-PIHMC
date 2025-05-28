#ifndef CPIHMC_RXN_COORD_DIFF_ASYM_H
#define CPIHMC_RXN_COORD_DIFF_ASYM_H

#include "rxn_coord_diff.h"

namespace cpihmc
{
    class rxn_coord_diff_asym : public rxn_coord_diff
    {
        private:
            const bool_t length_evolve() const final;
            const prec_t adjust_angle(const index_t) const;
            const bool_t angle_evolve() const final;
        public:
            explicit rxn_coord_diff_asym(const std::array<dof, 3> &RxnDOFs, const prec_t LengthWidth, const std::array<prec_t, 2> AngleWidths, const mc_mbr_pck &McMbrPck, const prec_t RxnCoordVal):rxn_coord_diff(RxnDOFs, LengthWidth, AngleWidths, McMbrPck, RxnCoordVal){}
            ~rxn_coord_diff_asym() = default;
            const prec_t get_mean_force_left() const final;
    };
}

#endif