#include "elec_num.h"
using namespace cpihmc;
using namespace std;

const bool_t cpihmc::elec_num::mc_evolve()
{
    bool_t Accept = true;
    const prec_t TempPotEng = Box->PotEng;
    const prec_t TempElecNum = Box->ElecNum;
    prec_t Prefactor = 1.0;
    pair<prec_t, prec_t> CurrRange, TrialRange;
    CurrRange.first = max(Box->ElecNum-0.5*ElecNumWidth, ElecNumRange.first);
    CurrRange.second = min(Box->ElecNum+0.5*ElecNumWidth, ElecNumRange.second);
    Box->ElecNum = CurrRange.first + Rand.uniform()*(CurrRange.second-CurrRange.first);
    TrialRange.first = max(Box->ElecNum-0.5*ElecNumWidth, ElecNumRange.first);
    TrialRange.second = min(Box->ElecNum+0.5*ElecNumWidth, ElecNumRange.second);
    Prefactor = (CurrRange.second-CurrRange.first) / (TrialRange.second-TrialRange.first);
    if (Box->ElecNum != TempElecNum)
    {
        Box->PotEng = Pot->infer(Box);
        if (Rand.uniform() > min(1.0, Prefactor*exp(((Box->ElecNum-TempElecNum)*Mu-(Box->PotEng-TempPotEng))/(k_B*Box->Temp))))
        {
            Box->ElecNum = TempElecNum;
            Temp->force_temp2current();
            Box->PotEng = TempPotEng;
            Accept = false;
        }
        else
        {
            Temp->force_current2temp();
            Box->update_total_energy();
        }
    }
    return Accept;
}