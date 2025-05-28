#ifndef CPIHMC_VIRT_ATOM_H
#define CPIHMC_VIRT_ATOM_H

#include "hmc.h"
#include "dof.h"
#include "const.h"

namespace cpihmc
{
    class virt_atom : public hmc
    {
        private:
            size_t AtomNum;
            std::vector<atom *> AllAtoms;
            std::vector<atom> TFAtoms;
        private:
            void update_virt_atom_coordinates();
            const vec3<prec_t> get_transformed_force(const index_t) const;
            void update_atom_coordinates() const;
        public:
            dof_atom Centroid;
        public:
            explicit virt_atom(const std::vector<atom *> &AllAtoms, const prec_t Temp, const prec_t TimeStep, const prec_t MassScal, rand &Rand);
            ~virt_atom() = default;
            void random_velocities() final;
            void hmc_evolve_forth_before() final;
            void hmc_evolve_forth_after() final;
            void hmc_evolve_back() final;
            void update_atom_velocities() const;
    };
}

#endif