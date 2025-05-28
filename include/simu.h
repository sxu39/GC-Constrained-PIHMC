#ifndef CPIHMC_SIMU_H
#define CPIHMC_SIMU_H

#include "pot.h"
#include "output.h"
#include "md.h"
#include "pimc.h"
#include "elec_num.h"
#include "rxn_coord.h"
#include "rigid_body.h"
#include "virt_atom.h"
#include "hard_bdry.h"
#include "model_devi.h"

namespace cpihmc
{
    class simu
    {
        private:
            const input &Input;
            output &Output;
            rand &Rand;
            box &Box;
            pot *Pot;
            md *MD;
            pimc *PIMC;
            elec_num *ElecNumEvolve;
            std::vector<index_t> HMCAtomIndexes;
            std::vector<rigid_body> RigidBodies;
            std::vector<virt_atom> VirtAtoms;
            std::vector<rxn_coord *> RxnCoords;
            temp *Temp;
            hard_bdry *HardBdry;
            std::array<prec_t, 3> Ratios;
            void (simu::*SingleEvolve)();
            model_devi *ModelDevi;
        private:
            static void adjust_rxn_coord(rxn_coord_info &, const index_t);
            const std::unordered_map<index_t, std::vector<index_t>> preset_rxn_coords(std::vector<rxn_coord_info> &) const;
            void set_up_rxn_coord(const mc_mbr_pck &);
            void current2temp();
            void temp2current();
            void random_velocities();
            void md_evolve();
            void hmc_evolve();
            rxn_coord * const random_rxn_coord() const;
            void mc_evolve(mc * const, void (temp::*CoordCurrent2Temp)()=nullptr);
            void chmc_evolve();
        public:
            explicit simu(const input &, output &, rand &, box &);
            ~simu();
            void evolve();
    };
}

#endif