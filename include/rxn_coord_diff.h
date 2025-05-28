#ifndef CPIHMC_RXN_COORD_DIFF_H
#define CPIHMC_RXN_COORD_DIFF_H

#include "rxn_coord.h"

namespace cpihmc
{
    enum {Length, Theta1, Phi1, Theta2, Phi2}; // The type of MC step
    class rxn_coord_diff : public rxn_coord
    {
        protected:
            index_t Type;
            std::array<dof, 3> RxnDOFs;
            prec_t LengthWidth;
            std::array<prec_t, 2> AngleWidths;
        private:
            virtual const bool_t length_evolve() const;
            const prec_t adjust_angle(const index_t) const;
            virtual const bool_t angle_evolve() const;
        public:
            explicit rxn_coord_diff(const std::array<dof, 3> &RxnDOFs, const prec_t LengthWidth, const std::array<prec_t, 2> AngleWidths, const mc_mbr_pck &McMbrPck, const prec_t RxnCoordVal):Type(Length), RxnDOFs(RxnDOFs), LengthWidth(LengthWidth), AngleWidths(AngleWidths), rxn_coord(McMbrPck, RxnCoordVal){}
            virtual ~rxn_coord_diff() = default;
            const bool_t mc_evolve() final;
            const prec_t get_rxn_coord_val() const final;
            const prec_t get_mean_force_left() const override;
            const prec_t get_mean_force_right() const final;
    };
}

#endif