#include "rxn_coord_diff.h"
using namespace cpihmc;
using namespace std;

const bool_t cpihmc::rxn_coord_diff::length_evolve() const
{
    bool_t Accept = false;
    const prec_t TempPotEng = Box->PotEng; // record current potential energy
    vec3<prec_t> LeftVec = RxnDOFs[0].Coord - RxnDOFs[1].Coord;
    vec3<prec_t> RightVec = RxnDOFs[2].Coord - RxnDOFs[1].Coord;
    prec_t CurrLeftDist = LeftVec.norm(), CurrRightDist = RightVec.norm();
    prec_t NewRightDist = CurrRightDist + (Rand.uniform()-0.5) * LengthWidth;
    prec_t NewLeftDist = NewRightDist + RxnCoordVal;
    RxnDOFs[0].Coord = RxnDOFs[1].Coord + NewLeftDist * LeftVec / CurrLeftDist;
    RxnDOFs[2].Coord = RxnDOFs[1].Coord + NewRightDist * RightVec / CurrRightDist;
    for (auto &Bead : Box->Beads) Bead.back_centroid();
    const prec_t Prefactor = pow(NewLeftDist*NewRightDist/CurrLeftDist/CurrRightDist, 2);
    if (HardBdry->judge_accept())
    {
        Box->PotEng = Pot->infer(Box);
        if (Rand.uniform() > min(1.0, Prefactor*exp(-(Box->PotEng-TempPotEng)/(k_B*Box->Temp))))
        {
            RxnDOFs[0].Coord = RxnDOFs[1].Coord + (CurrRightDist+RxnCoordVal) * LeftVec / CurrLeftDist;
            RxnDOFs[2].Coord = RxnDOFs[1].Coord + RightVec;
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
        RxnDOFs[0].Coord = RxnDOFs[1].Coord + (CurrRightDist+RxnCoordVal) * LeftVec / CurrLeftDist;
        RxnDOFs[2].Coord = RxnDOFs[1].Coord + RightVec;
        for (auto &Bead : Box->Beads) Bead.back_centroid();
    }
    return Accept;
}

const prec_t cpihmc::rxn_coord_diff::adjust_angle(const index_t TagDOFIndex) const
{
    vec3<prec_t> DistVec = RxnDOFs[TagDOFIndex].Coord - RxnDOFs[1].Coord;
    prec_t SinTheta[2];
    prec_t CurrTheta, NewTheta, CurrPhi, NewPhi;
    switch (Type)
    {
        case Theta1:
        case Theta2:
            SinTheta[0] = sqrt(1-pow(DistVec.z/DistVec.norm(), 2));
            CurrTheta = acos(DistVec.z/DistVec.norm());
            NewTheta = CurrTheta + (Rand.uniform()-0.5) * AngleWidths[0];
            SinTheta[1] = sin(NewTheta);
            RxnDOFs[TagDOFIndex].Coord = vec3<prec_t>{RxnDOFs[1].Coord.x + DistVec.x / SinTheta[0] * SinTheta[1], 
                                  RxnDOFs[1].Coord.y + DistVec.y / SinTheta[0] * SinTheta[1],
                                  RxnDOFs[1].Coord.z + DistVec.norm() * cos(NewTheta)};
            return SinTheta[1] / SinTheta[0];
            break;
        case Phi1:
        case Phi2:
            CurrTheta = acos(DistVec.z/DistVec.norm());
            CurrPhi = acos(DistVec.x/DistVec.norm()/sin(CurrTheta));
            NewPhi = CurrPhi + (Rand.uniform()-0.5) * AngleWidths[1];
            RxnDOFs[TagDOFIndex].Coord = vec3<prec_t>{RxnDOFs[1].Coord.x + DistVec.x * cos(NewPhi) / cos(CurrPhi), 
                                  RxnDOFs[1].Coord.y + DistVec.y * sin(NewPhi) / sin(CurrPhi),
                                  RxnDOFs[TagDOFIndex].Coord.z};
        default:
            return 1.0;
            break;
    }
}

const bool_t cpihmc::rxn_coord_diff::angle_evolve() const
{
    bool_t Accept = false;
    const prec_t TempPotEng = Box->PotEng; // record current potential energy
    const index_t TagDOFIndex = 2 * ((Type-1) / 2);
    const vec3<prec_t> TempCoord = RxnDOFs[TagDOFIndex].Coord;
    const prec_t Prefactor = adjust_angle(TagDOFIndex);
    for (auto &Bead : Box->Beads) Bead.back_centroid();
    if (HardBdry->judge_accept())
    {
        Box->PotEng = Pot->infer(Box);
        if (Rand.uniform() > min(1.0, Prefactor*exp(-(Box->PotEng-TempPotEng)/(k_B*Box->Temp))))
        {
            RxnDOFs[TagDOFIndex].Coord = TempCoord;
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
        RxnDOFs[TagDOFIndex].Coord = TempCoord;
        for (auto &Bead : Box->Beads) Bead.back_centroid();
    }
    return Accept;
}

const bool_t cpihmc::rxn_coord_diff::mc_evolve()
{
    Type = Rand.uniform_int(5);
    if (Type == Length) return length_evolve();
    else return angle_evolve();
}

const prec_t cpihmc::rxn_coord_diff::get_rxn_coord_val() const {return (RxnDOFs[0].Coord - RxnDOFs[1].Coord).norm() - (RxnDOFs[2].Coord - RxnDOFs[1].Coord).norm();}

const prec_t cpihmc::rxn_coord_diff::get_mean_force_left() const
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

const prec_t cpihmc::rxn_coord_diff::get_mean_force_right() const
{
    prec_t MeanForce = 0.0;

    // derivative of potential energy term
    vec3<prec_t> Direction = RxnDOFs[2].Coord - RxnDOFs[1].Coord;
    const prec_t LengthInv = 1.0 / Direction.norm();
    MeanForce += Direction * LengthInv * RxnDOFs[2].get_force();

    // Jacobian term
    MeanForce += 2.0 * k_B * Box->Temp * LengthInv;
    return MeanForce;
}