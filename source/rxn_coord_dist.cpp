#include "rxn_coord_dist.h"
using namespace cpihmc;
using namespace std;

const prec_t cpihmc::rxn_coord_dist::adjust_angle() const
{
    const vec3<prec_t> CurrDistVec = RxnDOFs[1].Coord - RxnDOFs[0].Coord;
    const prec_t CurrSinTheta = sqrt(1-pow(CurrDistVec.z/CurrDistVec.norm(), 2));

    vec3<prec_t> BallMove;
    prec_t Norm;
    do
    {
        BallMove.x = Rand.normal();
        BallMove.y = Rand.normal();
        BallMove.z = Rand.normal();
        Norm = BallMove.norm();
    } while (Norm == 0 || isinf(Norm)); // if the `BallMove` is 0 or has infinity number, produce `BallMove` again.

    // change the `BallMove` length
    const prec_t Length = Rand.uniform() * Radius;
    BallMove *= Length / Norm;

    RxnDOFs[1].Coord += BallMove;
    const vec3<prec_t> NewDistVec = RxnDOFs[1].Coord - RxnDOFs[0].Coord;
    const prec_t NewSinTheta = sqrt(1-pow(NewDistVec.z/NewDistVec.norm(), 2));
    RxnDOFs[1].Coord = RxnDOFs[0].Coord + NewDistVec * RxnCoordVal / NewDistVec.norm();
    return NewSinTheta / CurrSinTheta;
}

const bool_t cpihmc::rxn_coord_dist::mc_evolve()
{
    bool_t Accept = false;
    const prec_t TempPotEng = Box->PotEng;
    const vec3<prec_t> TempCoord = RxnDOFs[1].Coord;
    const prec_t Prefactor = adjust_angle();
    for (auto &Bead : Box->Beads) Bead.back_centroid();
    if (HardBdry->judge_accept())
    {
        Box->PotEng = Pot->infer(Box);
        if (Rand.uniform() > min(1.0, Prefactor*exp(-(Box->PotEng-TempPotEng)/(k_B*Box->Temp))))
        {
            RxnDOFs[1].Coord = TempCoord;
            for (auto &Bead : Box->Beads) Bead.back_centroid();
            Temp->force_temp2current();
            Box->PotEng = TempPotEng;
        }
        else
        {
            Temp->force_current2temp();
            Box->update_total_energy();
            Accept = true;
        }
    }
    else
    {
        RxnDOFs[1].Coord = TempCoord;
        for (auto &Bead : Box->Beads) Bead.back_centroid();
    }
    return Accept;
}

const prec_t cpihmc::rxn_coord_dist::get_rxn_coord_val() const {return (RxnDOFs[1].Coord - RxnDOFs[0].Coord).norm();}

const prec_t cpihmc::rxn_coord_dist::get_mean_force_left() const
{
    prec_t MeanForce = 0.0;

    // derivative of potential energy term
    vec3<prec_t> Direction = RxnDOFs[0].Coord - RxnDOFs[1].Coord;
    const prec_t LengthInv = 1.0 / Direction.norm();
    MeanForce -= Direction * LengthInv * RxnDOFs[0].get_force();

    // Jacobian term
    MeanForce -= 2.0 * k_B * Box->Temp * LengthInv;
    return MeanForce;
}

const prec_t cpihmc::rxn_coord_dist::get_mean_force_right() const
{
    prec_t MeanForce = 0.0;

    // derivative of potential energy term
    vec3<prec_t> Direction = RxnDOFs[1].Coord - RxnDOFs[0].Coord;
    const prec_t LengthInv = 1.0 / Direction.norm();
    MeanForce -= Direction * LengthInv * RxnDOFs[1].get_force();

    // Jacobian term
    MeanForce -= 2.0 * k_B * Box->Temp * LengthInv;
    return MeanForce;
}