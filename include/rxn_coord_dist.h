#ifndef CPIHMC_RXN_COORD_DIST_H
#define CPIHMC_RXN_COORD_DIST_H

#include "rxn_coord.h"

namespace cpihmc
{
    class rxn_coord_dist : public rxn_coord
    {
        private:
            std::array<dof, 2> RxnDOFs;
            prec_t Radius;
        private:
            const prec_t adjust_angle() const;
        public:
            explicit rxn_coord_dist(const std::array<dof, 2> &RxnDOFs, const prec_t Radius, const mc_mbr_pck &McMbrPck, const prec_t RxnCoordVal):RxnDOFs(RxnDOFs), Radius(Radius), rxn_coord(McMbrPck, RxnCoordVal){}
            ~rxn_coord_dist() = default;
            const bool_t mc_evolve() final;
            const prec_t get_rxn_coord_val() const final;
            const prec_t get_mean_force_left() const final;
            const prec_t get_mean_force_right() const final;
    };
}

#endif