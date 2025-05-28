#include "pimc.h"
using namespace cpihmc;
using namespace std;

void cpihmc::pimc::u2x(bead * const Bead, const index_t ChangeBeadIndex) const
{
    const size_t P = Box->NBead;
    for (index_t i = ChangeBeadNum ; i > 0 ; --i)
    {
        vec3<prec_t> u;
        u.x = sqrt(i*pow(hbar, 2)/((i+1)*Bead->Mass*P*k_B*Box->Temp)) * Rand.normal();
        u.y = sqrt(i*pow(hbar, 2)/((i+1)*Bead->Mass*P*k_B*Box->Temp)) * Rand.normal();
        u.z = sqrt(i*pow(hbar, 2)/((i+1)*Bead->Mass*P*k_B*Box->Temp)) * Rand.normal();
        Bead->Coords[(ChangeBeadIndex+i)%P] = u + ((prec_t)i*Bead->Coords[(ChangeBeadIndex+i+1)%P]+Bead->Coords[ChangeBeadIndex%P]) * (1/(prec_t)(i+1));
    }
}

const bool_t cpihmc::pimc::mc_evolve()
{
    bool_t Accept = false;
    prec_t TempPotEng = Box->PotEng;
    prec_t TempQtmKinEng = Box->QtmKinEng;
    const size_t BeadNum = Box->Beads.size();
    for (index_t i = 0 ; i < BeadNum ; ++i) TempBeadCoords[i] = Box->Beads[i].Coords;
    
    index_t ChangeBeadIndex = Rand.uniform_int(Box->NBead);
    for (index_t i = 0 ; i < BeadNum ; ++i)
    {
        u2x(&(Box->Beads[i]), ChangeBeadIndex);
        Box->Beads[i].back_centroid();
    }
    if (HardBdry->judge_accept())
    {
        Box->PotEng = Pot->infer(Box);
        if (Rand.uniform() > min(1.0, exp(-(Box->PotEng-TempPotEng)/(k_B*Box->Temp))))
        {
            for (index_t i = 0 ; i < BeadNum ; ++i) Box->Beads[i].Coords = TempBeadCoords[i];
            Temp->force_temp2current();
            Box->PotEng = TempPotEng;
            Box->QtmKinEng = TempQtmKinEng;
        }
        else
        {
            Temp->force_current2temp();
            Box->update_total_energy();
            Accept = true;
        }
    }
    else for (index_t i = 0 ; i < BeadNum ; ++i) Box->Beads[i].Coords = TempBeadCoords[i];
    return Accept;
}