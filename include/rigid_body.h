#ifndef CPIHMC_RIGID_BODY_H
#define CPIHMC_RIGID_BODY_H

#include "hmc.h"
#include "dof.h"
#include "const.h"

namespace cpihmc
{
    class rigid_body : public hmc
    {
        private:
            std::vector<dof> RigidBodyDOFs;
            const prec_t TotalMass;
        private:
            inline const prec_t get_total_mass() const;
            inline const vec3<prec_t> get_total_force() const;
        public:
            explicit rigid_body(const std::vector<dof> &RigidBodyDOFs, const prec_t Temp, const prec_t TimeStep, const prec_t MassScal, rand &Rand):RigidBodyDOFs(RigidBodyDOFs), TotalMass(get_total_mass()), hmc(Temp, TimeStep, MassScal, Rand){}
            ~rigid_body() = default;
            void random_velocities() final;
            void hmc_evolve_forth_before() final;
            void hmc_evolve_forth_after() final;
            void hmc_evolve_back() final;
    };
}

const cpihmc::prec_t cpihmc::rigid_body::get_total_mass() const
{
    prec_t TotalMass = 0.0;
    for (const auto &DOF : RigidBodyDOFs) TotalMass += DOF.Mass * DOF.AtomNum; // the mass of the DOF should be multiplied by the number of atom(s) involved in the DOF
    return TotalMass;
}

const cpihmc::vec3<cpihmc::prec_t> cpihmc::rigid_body::get_total_force() const
{
    vec3<prec_t> TotalForce;
    for (const auto &DOF : RigidBodyDOFs) TotalForce += DOF.get_force();
    return TotalForce;
}

#endif