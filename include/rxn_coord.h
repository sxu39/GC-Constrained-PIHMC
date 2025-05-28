#ifndef CPIHMC_RXN_COORD_H
#define CPIHMC_RXN_COORD_H

#include "mc.h"
#include "dof.h"

namespace cpihmc
{
    class rxn_coord : public mc
    {
        protected:
            prec_t RxnCoordVal;
        public:
            explicit rxn_coord(const mc_mbr_pck &McMbrPck, const prec_t RxnCoordVal):RxnCoordVal(RxnCoordVal), mc(McMbrPck){}
            virtual ~rxn_coord() = default;
            void set_rxn_coord_val(){RxnCoordVal = get_rxn_coord_val();}
            virtual const prec_t get_rxn_coord_val() const = 0;
            virtual const prec_t get_mean_force_left() const = 0; // calculate the mean force of left atom
            virtual const prec_t get_mean_force_right() const = 0; // calculate the mean force of right atom
    };
}

#endif