#include "rigid_body.h"
using namespace cpihmc;

void cpihmc::rigid_body::random_velocities()
{
    vec3<prec_t> Vel;
    Vel.x = sqrt(Temp*k_B/(TotalMass*MassScal)) * Rand.normal();
    Vel.y = sqrt(Temp*k_B/(TotalMass*MassScal)) * Rand.normal();
    Vel.z = sqrt(Temp*k_B/(TotalMass*MassScal)) * Rand.normal();
    for (auto &DOF : RigidBodyDOFs) DOF.Vel = Vel;
}

void cpihmc::rigid_body::hmc_evolve_forth_before()
{
    const vec3<prec_t> TotalForce = get_total_force();
    for (auto &DOF : RigidBodyDOFs)
    {
        DOF.Vel += TotalForce * TimeStep * 0.5 / (TotalMass * MassScal);
        DOF.Coord += DOF.Vel * TimeStep;
    }
}

void cpihmc::rigid_body::hmc_evolve_forth_after()
{
    const vec3<prec_t> TotalForce = get_total_force();
    for (auto &DOF : RigidBodyDOFs) DOF.Vel += TotalForce * TimeStep * 0.5 / (TotalMass * MassScal);
}

void cpihmc::rigid_body::hmc_evolve_back(){for (auto &DOF : RigidBodyDOFs) DOF.synchronize();}