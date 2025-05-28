#ifndef CPIHMC_TEMP_H
#define CPIHMC_TEMP_H

#include "box.h"

namespace cpihmc
{
    class temp
    {
        private:
            box * const Box;
            std::vector<vec3<prec_t>> TempCoords;
            std::vector<vec3<prec_t>> TempForces;
            std::vector<index_t> RxnCoordIndexes;
        public:
            explicit temp(box * const Box):Box(Box), TempCoords(Box->NAtoms), TempForces(Box->NAtoms){coord_current2temp(); force_current2temp();}
            ~temp() = default;
            void set_rxn_coord_indexes(const std::vector<index_t> &RxnCoordIndexes){this->RxnCoordIndexes = RxnCoordIndexes;}
            void coord_current2temp(){for (index_t i = 0 ; i < Box->NAtoms ; ++i) TempCoords[i] = Box->Atoms[i].Coord;}
            void coord_temp2current() const
            {
                for (index_t i = 0 ; i < Box->NAtoms ; ++i) Box->Atoms[i].Coord = TempCoords[i];
                for (auto &Bead : Box->Beads) Bead.back_centroid();
            }
            void coord_current2temp_rxn_coord(){for (const auto &i : RxnCoordIndexes) TempCoords[i] = Box->Atoms[i].Coord;}
            void force_current2temp(){for (index_t i = 0 ; i < Box->NAtoms ; ++i) TempForces[i] = Box->Atoms[i].Force;}
            void force_temp2current() const {for (index_t i = 0 ; i < Box->NAtoms ; ++i) Box->Atoms[i].Force = TempForces[i];}
    };
}

#endif