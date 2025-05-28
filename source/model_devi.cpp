#include "model_devi.h"
using namespace cpihmc;
using namespace std;

cpihmc::model_devi::model_devi(box * const Box, const vector<string> &ModelNames, const string &FileName):Box(Box), ModelNum(ModelNames.size()), OutFile(FileName)
{
    Models.resize(ModelNum);
    for (index_t i = 0 ; i < ModelNum ; ++i) Models[i].init(ModelNames[i]);
    OutFile << scientific;
    OutFile << "#" << setw(12 - 1) << "step" << setw(18 + 1) << "max_devi_v" << setw(18 + 1) << "min_devi_v" << setw(18 + 1) << "avg_devi_v"
                                             << setw(18 + 1) << "max_devi_f" << setw(18 + 1) << "min_devi_f" << setw(18 + 1) << "avg_devi_f" << endl;
}

model_devi::~model_devi(){if (OutFile.is_open()) OutFile.close();}

void cpihmc::model_devi::compute_force_virials(vector<vector<prec_t>> &AllForces, vector<vector<prec_t>> &AllVirials, const vector<prec_t> &Coords, const vector<index_t> &AtomTypes, const vector<prec_t> &Cells)
{
    vector<prec_t> Energies;
    AllForces.resize(ModelNum);
    AllVirials.resize(ModelNum);
    for (index_t i = 0 ; i < ModelNum ; ++i)
    {
        if (Models[i].dim_fparam()) Models[i].compute(Energies, AllForces[i], AllVirials[i], Coords, AtomTypes, Cells, vector<prec_t>(Box->NBead, Box->ElecNum));
        else Models[i].compute(Energies, AllForces[i], AllVirials[i], Coords, AtomTypes, Cells);
    }
}

template <class T>
void cpihmc::model_devi::compute_avg(vector<T> &TempAvg, const vector<vector<T>> &Vec) const
{
    assert(Vec.size() == ModelNum);
    if (ModelNum == 0) return;
    TempAvg.resize(Vec[0].size(), (T)0.0);
    for (index_t ii = 0 ; ii < ModelNum ; ++ii) for (index_t jj = 0 ; jj < TempAvg.size() ; ++jj) TempAvg[jj] += Vec[ii][jj];
    for (auto &Val : TempAvg) Val /= (T)ModelNum;
}

template void cpihmc::model_devi::compute_avg<prec_t>(vector<prec_t> &, const vector<vector<prec_t>> &) const;

void cpihmc::model_devi::ana_st(prec_t &Max, prec_t &Min, prec_t &Avg, const vector<prec_t> &Vec)
{
    if (Vec.size())
    {
        Max = numeric_limits<prec_t>::min();
        Min = numeric_limits<prec_t>::max();
        Avg = 0.0;
        for (const auto &Val : Vec)
        {
            if (Val > Max) Max = Val;
            if (Val < Min) Min = Val;
            Avg += Val;
        }
        Avg /= Vec.size();
    }
}

template <class T>
void cpihmc::model_devi::compute_std(vector<T> &Std, const vector<T> &Avg, const vector<vector<T>>& Vec, const index_t Stride) const
{
    assert(Vec.size() == ModelNum);
    if (ModelNum == 0) return;

    size_t DOFNum = Avg.size();
    size_t LocNum = DOFNum / Stride;
    assert(LocNum * Stride == DOFNum);

    Std.resize(LocNum, (T)0.0);

    for (index_t ii = 0 ; ii < ModelNum ; ++ii)
        for (index_t jj = 0 ; jj < LocNum ; ++jj)
        {
            const T *TempVal = &(Vec[ii][jj*Stride]);
            const T *TempAvg = &(Avg[jj*Stride]);
            for (index_t dd = 0 ; dd < Stride ; ++dd)
            {
                const T Diff = TempVal[dd] - TempAvg[dd];
                Std[jj] += Diff * Diff;
            }
        }

    for (index_t jj = 0 ; jj < LocNum ; ++jj) Std[jj] = sqrt(Std[jj] / (T)ModelNum);
}

template void cpihmc::model_devi::compute_std<prec_t>(vector<prec_t> &, const vector<prec_t> &, const vector<vector<prec_t>> &, const index_t) const;

void cpihmc::model_devi::cope_with_force(prec_t &ForceMax, prec_t &ForceMin, prec_t &ForceAvg, const vector<vector<prec_t>> &AllForces) const
{
    vector<prec_t> TempForceAvg;
    vector<prec_t> ForceStd;
    compute_avg(TempForceAvg, AllForces);
    compute_std(ForceStd, TempForceAvg, AllForces, 3);
    ana_st(ForceMax, ForceMin, ForceAvg, ForceStd);
    ForceMax *= eV_d_Ang2Force;
    ForceMin *= eV_d_Ang2Force;
    ForceAvg *= eV_d_Ang2Force;
}

void cpihmc::model_devi::cope_with_virial(prec_t &VirialMax, prec_t &VirialMin, prec_t &VirialAvg, const vector<vector<prec_t>> &AllVirials) const
{
    vector<prec_t> TempVirialAvg;
    vector<prec_t> VirialStd;
    compute_avg(TempVirialAvg, AllVirials);
    compute_std(VirialStd, TempVirialAvg, AllVirials, 1);
    for (const auto &Val : VirialStd)
    {
        if (Val > VirialMax) VirialMax = Val;
        if (Val < VirialMin) VirialMin = Val;
        VirialAvg += Val * Val;
    }
    VirialAvg = sqrt(VirialAvg / (9*Box->NBead));
    VirialMax *= eV2Energy;
    VirialMin *= eV2Energy;
    VirialAvg *= eV2Energy;
}

void cpihmc::model_devi::compute_model_devi(const index_t StepIndex)
{
    // collect all beads structures in the format required for DP inference
    const index_t RealAtomNum = Box->AtomIndexes.size();
    vector<prec_t> Coords(Box->NBead*RealAtomNum*3);
    vector<prec_t> Cells;
    vector<index_t> AtomTypes(RealAtomNum);
    deep_pot::provide_structures(Coords, AtomTypes, Cells, Box);

    vector<vector<prec_t>> AllForces;
    vector<vector<prec_t>> AllVirials;
    compute_force_virials(AllForces, AllVirials, Coords, AtomTypes, Cells);

    prec_t ForceMin = numeric_limits<prec_t>::max(), ForceMax = 0.0, ForceAvg = 0.0;
    cope_with_force(ForceMax, ForceMin, ForceAvg, AllForces);

    for (index_t kk = 0 ; kk < ModelNum ; ++kk) for (index_t ii = 0 ; ii < 9 * Box->NBead ; ++ii) AllVirials[kk][ii] /= RealAtomNum;
    prec_t VirialMin = numeric_limits<prec_t>::max(), VirialMax = 0.0, VirialAvg = 0.0;
    cope_with_virial(VirialMax, VirialMin, VirialAvg, AllVirials);

    OutFile << setw(12) << StepIndex << " " << setw(18) << VirialMax << " " << setw(18) << VirialMin << " " << setw(18) << VirialAvg
                                     << " " << setw(18) << ForceMax << " " << setw(18) << ForceMin << " " << setw(18) << ForceAvg << endl;
}