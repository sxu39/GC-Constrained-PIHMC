#include "simu.h"
#include "ml_pot.h"
#include "md_nve.h"
#include "rxn_coord_dist.h"
#include "rxn_coord_diff.h"
#include "rxn_coord_diff_asym.h"
#include "wall.h"
#include <algorithm>
using namespace cpihmc;
using namespace std;

cpihmc::simu::simu(const input &Input, output &Output, rand &Rand, box &Box):Input(Input), Output(Output), Rand(Rand), Box(Box), Pot(nullptr), MD(nullptr), PIMC(nullptr), ElecNumEvolve(nullptr), Temp(nullptr), HardBdry(nullptr), Ratios({0.0, 0.0, 0.0}), SingleEvolve(nullptr), ModelDevi(nullptr)
{
    if (Input.get_pot_type() == "DP") Pot = new deep_pot{Input.get_deep_pot_model()};
    this->Box.PotEng = Pot->infer(&Box);
    this->Box.update_total_energy();
    if (Input.get_simu_type() == "HMC" || Input.get_simu_type() == "CHMC")
    {
        Temp = new temp{&Box};
        HardBdry = new wall{&Box, Input.get_wall()};
        HMCAtomIndexes = Box.AtomIndexes;
        if (Input.get_simu_type() == "CHMC")
        {
            const mc_mbr_pck McMbrPck{&Box, Pot, Temp, HardBdry, Rand};
            set_up_rxn_coord(McMbrPck);
            if (Input.get_elec_num_ratio())
            {
                if (Box.ElecNum < Input.get_elec_num_range().first) throw invalid_argument("Current electron number "+to_string(Box.ElecNum)+" is lower than the lower bound.");
                else if (Box.ElecNum > Input.get_elec_num_range().second) throw invalid_argument("Current electron number "+to_string(Box.ElecNum)+" is higher than the upper bound.");
                ElecNumEvolve = new elec_num{McMbrPck, Input.get_mu(), Input.get_elec_num_range(), Input.get_elec_num_width()};
                Ratios[0] = Input.get_elec_num_ratio();
                Ratios[1] = Ratios[0];
            }
            if (Box.NBead > 1)
            {
                PIMC = new pimc{McMbrPck, Input.get_n_change_bead()};
                const prec_t RepeatTimes = (prec_t)Input.get_n_bead() / (prec_t)Input.get_n_change_bead();
                Ratios[1] = Ratios[0] + RepeatTimes / (1+RepeatTimes) * (1-Ratios[0]);
            }
            Ratios[2] = 1 - Input.get_hybrid_monte_carlo_ratio() * (1-Ratios[1]);
            SingleEvolve = &simu::chmc_evolve;
        }
        else SingleEvolve = &simu::hmc_evolve;
    }
    else if (Input.get_simu_type() == "MD")
    {
        if (abs(Box.MassScal-1.0) > 1e-9) throw invalid_argument("The mass scaling coefficient cannot be used in the MD simulation.");
        if (Box.NBead > 1) throw invalid_argument("The path integral formalism cannot be used in the MD simulation.");
        if (Input.get_wall().size()) throw invalid_argument("The hard wall cannot be used in the MD simulation.");
        MD = new md_nve{&this->Box, Pot, Input.get_time_step()};
        SingleEvolve = &simu::md_evolve;
        random_velocities();
        Box.update_total_energy(true);
    }
    if (Input.get_model_devi_deep_pot_models().size()) ModelDevi = new model_devi{&this->Box, Input.get_model_devi_deep_pot_models(), Input.get_model_devi_file()};
}

cpihmc::simu::~simu()
{
    if (Pot) delete Pot;
    if (MD) delete MD;
    for (auto &RxnCoord : RxnCoords) delete RxnCoord;
    if (PIMC) delete PIMC;
    if (ElecNumEvolve) delete ElecNumEvolve;
    if (Temp) delete Temp;
    if (HardBdry) delete HardBdry;
    if (ModelDevi) delete ModelDevi;
}

void cpihmc::simu::adjust_rxn_coord(rxn_coord_info &Info, const index_t Index)
{
    if (Info.type == "DIST"){if (Index == 1) swap(Info.atom_indexes[0], Info.atom_indexes[1]);}
    else if (Info.type == "DIFF") if (Index == 2){Info.type = "DIFF_ASYM"; swap(Info.atom_indexes[0], Info.atom_indexes[2]);}
}

const unordered_map<index_t, vector<index_t>> cpihmc::simu::preset_rxn_coords(vector<rxn_coord_info> &RxnCoordInfo) const
{
    const size_t RxnCoordNum = RxnCoordInfo.size();
    vector<bool_t> Intsct(RxnCoordNum, false); // to save whether each reaction coordinate intersects with another reaction coordinate
    vector<index_t> IntsctIndexes(RxnCoordNum);
    for (index_t i = 0 ; i < RxnCoordNum ; ++i)
        for (index_t j = i + 1 ; j < RxnCoordNum ; ++j)
            for (index_t I = 0 ; I < RxnCoordInfo[i].atom_indexes.size() ; ++I)
                for (index_t J = 0 ; J < RxnCoordInfo[j].atom_indexes.size() ; ++J)
                    if (RxnCoordInfo[i].atom_indexes[I] == RxnCoordInfo[j].atom_indexes[J])
                    {
                        if (Intsct[i]){if (RxnCoordInfo[i].atom_indexes[I] != IntsctIndexes[i]) throw invalid_argument("A reaction coordinate cannot intersect with other reaction coordinate with more than one atom.");}
                        else
                        {
                            Intsct[i] = true;
                            IntsctIndexes[i] = RxnCoordInfo[i].atom_indexes[I];
                            adjust_rxn_coord(RxnCoordInfo[i], I);
                        }
                        if (Intsct[j]){if (RxnCoordInfo[j].atom_indexes[J] != IntsctIndexes[j]) throw invalid_argument("A reaction coordinate cannot intersect with other reaction coordinate with more than one atom.");}
                        else
                        {
                            Intsct[j] = true;
                            IntsctIndexes[j] = RxnCoordInfo[j].atom_indexes[J];
                            adjust_rxn_coord(RxnCoordInfo[j], J);
                        }
                    }
    unordered_map<index_t, vector<index_t>> RigidBodyIndexes;
    for (index_t i = 0 ; i < RxnCoordNum ; ++i)
    {
        if (Intsct[i])
        {
            if (RigidBodyIndexes.find(IntsctIndexes[i]) != RigidBodyIndexes.end())
                for (const auto &Index : RxnCoordInfo[i].atom_indexes){if (Index != IntsctIndexes[i]) RigidBodyIndexes[IntsctIndexes[i]].push_back(Index);}
            else
                RigidBodyIndexes.insert({IntsctIndexes[i], RxnCoordInfo[i].atom_indexes});
        }
        else RigidBodyIndexes.insert({RxnCoordInfo[i].atom_indexes[0], RxnCoordInfo[i].atom_indexes});
    }
    return RigidBodyIndexes;
}

void cpihmc::simu::set_up_rxn_coord(const mc_mbr_pck &McMbrPck)
{
    vector<index_t> RxnCoordIndexes;
    const vector<vector<index_t>> &VirtAtomIndexes = Input.get_virt_atom_indexes();
    for (const auto &VirtAtomIndex : VirtAtomIndexes)
    {
        vector<atom *> AllAtoms;
        for (const auto &Index : VirtAtomIndex) AllAtoms.push_back(&Box.Atoms[Index]);
        VirtAtoms.push_back(virt_atom{AllAtoms, Box.Temp, Input.get_time_step(), Box.MassScal, Rand});
        RxnCoordIndexes.insert(RxnCoordIndexes.end(), VirtAtomIndex.begin(), VirtAtomIndex.end());
    }
    vector<rxn_coord_info> RxnCoordInfo = Input.get_rxn_coord();
    const auto RigidBodyIndexes = preset_rxn_coords(RxnCoordInfo);
    for (const auto &Info : RxnCoordInfo)
    {
        vector<dof> RxnDOFs;
        for (const auto &Index : Info.atom_indexes)
        {
            if (Index < 0) RxnDOFs.push_back(dof{VirtAtoms[-1-Index].Centroid});
            else RxnDOFs.push_back(dof{Box.Atoms[Index]});
        }
        rxn_coord *RxnCoord = nullptr;
        if (Info.type == "DIST") RxnCoord = new rxn_coord_dist{{RxnDOFs[0], RxnDOFs[1]}, Info.params[0], McMbrPck, Info.value};
        if (Info.type == "DIFF") RxnCoord = new rxn_coord_diff{{RxnDOFs[0], RxnDOFs[1], RxnDOFs[2]}, Info.params[0], {Info.params[1], Info.params[2]}, McMbrPck, Info.value};
        if (Info.type == "DIFF_ASYM") RxnCoord = new rxn_coord_diff_asym{{RxnDOFs[0], RxnDOFs[1], RxnDOFs[2]}, Info.params[0], {Info.params[1], Info.params[2]}, McMbrPck, Info.value};
        if (RxnCoord){if (!Info.set_rxn_coord) RxnCoord->set_rxn_coord_val(); RxnCoords.push_back(RxnCoord);}
    }
    for (auto Iter = RigidBodyIndexes.begin() ; Iter != RigidBodyIndexes.end() ; ++Iter)
    {
        vector<dof> DOFs;
        for (const auto &Index : (*Iter).second)
        {
            if (Index < 0) DOFs.push_back(dof{VirtAtoms[-1-Index].Centroid});
            else
            {
                DOFs.push_back(dof{Box.Atoms[Index]});
                RxnCoordIndexes.push_back(Index);
            }
        }
        RigidBodies.push_back(rigid_body{DOFs, Box.Temp, Input.get_time_step(), Box.MassScal, Rand});
    }
    for (const auto &Index : RxnCoordIndexes) if (find(Box.Masks.begin(), Box.Masks.end(), Index) != Box.Masks.end()) throw invalid_argument("The removed atom with the index "+to_string(Index)+" cannot be involved in the constraint.");
    Temp->set_rxn_coord_indexes(RxnCoordIndexes);
    sort(RxnCoordIndexes.rbegin(), RxnCoordIndexes.rend());
    for (const auto &Index : RxnCoordIndexes) HMCAtomIndexes.erase(HMCAtomIndexes.begin()+Index);
}

void cpihmc::simu::current2temp(){Temp->coord_current2temp(); Temp->force_current2temp();}

void cpihmc::simu::temp2current(){Temp->coord_temp2current(); Temp->force_temp2current();}

void cpihmc::simu::random_velocities()
{
    for (const auto &i : Box.AtomIndexes)
    {
        if (Box.Atoms[i].Move.x) Box.Atoms[i].Vel.x = sqrt(Box.Temp*k_B/(Box.Atoms[i].Mass*Box.MassScal)) * Rand.normal();
        if (Box.Atoms[i].Move.y) Box.Atoms[i].Vel.y = sqrt(Box.Temp*k_B/(Box.Atoms[i].Mass*Box.MassScal)) * Rand.normal();
        if (Box.Atoms[i].Move.z) Box.Atoms[i].Vel.z = sqrt(Box.Temp*k_B/(Box.Atoms[i].Mass*Box.MassScal)) * Rand.normal();
    }
}

void cpihmc::simu::md_evolve(){MD->evolve();}

void cpihmc::simu::hmc_evolve()
{
    const prec_t TimeStep = Input.get_time_step();
    const size_t NEvolStep = Input.get_n_evol_step();
    random_velocities();
    for (auto &RigidBody : RigidBodies) RigidBody.random_velocities();
    for (auto &VirtAtoms : VirtAtoms) VirtAtoms.random_velocities();
    for (auto &VirtAtoms : VirtAtoms) VirtAtoms.update_atom_velocities();
    Box.update_total_energy(true);
    const prec_t TempKinEng = Box.KinEng;
    const prec_t TempPotEng = Box.PotEng;
    const prec_t TempTotEng = Box.TotEng;
    for (index_t step = 0 ; step < NEvolStep; ++step)
    {
        for (const auto &i : HMCAtomIndexes)
            if (Box.Atoms[i].Move.x * Box.Atoms[i].Move.y * Box.Atoms[i].Move.z)
            {
                Box.Atoms[i].Vel += Box.Atoms[i].Force * TimeStep * 0.5 / (Box.Atoms[i].Mass * Box.MassScal);
                Box.Atoms[i].Coord += Box.Atoms[i].Vel * TimeStep;
            }
        for (auto &RigidBody : RigidBodies) RigidBody.hmc_evolve_forth_before();
        for (auto &VirtAtom : VirtAtoms) VirtAtom.hmc_evolve_forth_before();
        for (auto &Bead : Box.Beads) Bead.back_centroid();
        Box.PotEng = Pot->infer(&Box);
        for (const auto &i : HMCAtomIndexes)
            if (Box.Atoms[i].Move.x * Box.Atoms[i].Move.y * Box.Atoms[i].Move.z)
                Box.Atoms[i].Vel += Box.Atoms[i].Force * TimeStep * 0.5 / (Box.Atoms[i].Mass * Box.MassScal);
        for (auto &RigidBody : RigidBodies) RigidBody.hmc_evolve_forth_after();
        for (auto &VirtAtom : VirtAtoms) VirtAtom.hmc_evolve_forth_after();
    }
    for (auto &VirtAtom : VirtAtoms) VirtAtom.update_atom_velocities();
    Box.update_total_energy(true);
    if (Rand.uniform() > min(1.0, exp(-(Box.TotEng-TempTotEng)/k_B/Box.Temp)) || !(HardBdry->judge_accept()))
    {
        temp2current();
        for (auto &RigidBody : RigidBodies) RigidBody.hmc_evolve_back();
        for (auto &VirtAtom : VirtAtoms) VirtAtom.hmc_evolve_back();
        Box.KinEng = TempKinEng;
        Box.PotEng = TempPotEng;
        Box.TotEng = TempTotEng;
        ++Output.RejectTimes;
    }
    else current2temp();
}

void cpihmc::simu::mc_evolve(mc * const MC, void (temp::*CoordCurrent2Temp)())
{
    const bool_t Accept = MC->mc_evolve();
    if (Accept){if (CoordCurrent2Temp) (Temp->*CoordCurrent2Temp)();}
    else ++Output.RejectTimes;
}

rxn_coord * const cpihmc::simu::random_rxn_coord() const {return RxnCoords[(RxnCoords.size() == 1) ? 0 : Rand.uniform_int(RxnCoords.size())];}

void cpihmc::simu::chmc_evolve()
{
    const prec_t RandNum = Rand.uniform();
    if (RandNum < Ratios[0]) mc_evolve(ElecNumEvolve);
    else if (RandNum < Ratios[1]) mc_evolve(PIMC);
    else if (RandNum < Ratios[2]) mc_evolve(random_rxn_coord(), &temp::coord_current2temp_rxn_coord);
    else hmc_evolve();
}

void cpihmc::simu::evolve()
{
    Output.BeginTime = time(0);
    Output.RejectTimes = 0;
    Output.update_physical_quantities(0, RxnCoords);
    Output.print_physical_quantities(0);
    Output.save_structure(0);
    if (ModelDevi) ModelDevi->compute_model_devi(0);
    size_t steps = Input.get_steps(), phy_quant_intvl = Input.get_phy_quant_intvl(), stru_intvl = Input.get_stru_intvl(), model_devi_intvl = Input.get_model_devi_intvl();
    for (index_t i = 0 ; i < steps ; ++i)
    {
        (this->*SingleEvolve)();
        Output.update_physical_quantities(i+1, RxnCoords);
        if ((i+1) % phy_quant_intvl == 0)
            Output.print_physical_quantities(i+1);
        if ((i+1) % stru_intvl == 0)
            Output.save_structure(i+1);
        if (ModelDevi)
            if ((i+1) % model_devi_intvl == 0)
                ModelDevi->compute_model_devi(i+1);
    }
    Output.EndTime = time(0);
    Output.simulation_summary();
}