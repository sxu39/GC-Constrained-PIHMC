#include "output.h"
using namespace cpihmc;
using namespace std;

unordered_map<string, string> cpihmc::output::PhyQuantName = {{"ke", "KinEng"}, {"pe", "PotEng"}, {"etotal", "TotEng"}, {"eint", "IntEng"}, 
                                                              {"pe_ave", "PotEngAve"}, {"pe_var", "PotEngVar"}, {"ne", "ElecNum"}};

cpihmc::output::output(const input &Input, const box &Box):Input(Input), Box(Box), OutFile(Input.get_phy_quant_file()), 
                                                           OutQuant({{"ke", 0.0}, {"pe", 0.0}, {"etotal", 0.0}, {"eint", 0.0}, 
                                                           {"pe_ave", 0.0}, {"pe_var", 0.0}, {"ne", 0.0}}), 
                                                           StepColWidth(max((size_t)floor(log10(Input.get_steps()))+1, (size_t)5)),
                                                           ColWidth(Input.get_col_width()), RejectTimes(0)
{
    const vector<rxn_coord_info> &RxnCoordInfo = Input.get_rxn_coord();
    if (RxnCoordInfo.size() == 1) insert_rxn_coord_quant();
    else for (index_t i = 0 ; i < RxnCoordInfo.size() ; ++i) insert_rxn_coord_quant("_"+to_string(i));
    OutFile.setf(ios::fixed);
    OutFile.setf(ios::showpoint);
    system("mkdir ALL_STRU");
    const vector<string> &OutPhyQuant = Input.get_out_phy_quant();
    OutFile << setw(StepColWidth) << "Steps";
    for (vector<string>::const_iterator it = OutPhyQuant.begin() ; it < OutPhyQuant.end() ; ++it)
        OutFile << setw(ColWidth) << PhyQuantName[*it];
    OutFile << "\n";
}

void cpihmc::output::insert_rxn_coord_quant(const std::string &Suffix)
{
    PhyQuantName.insert({"rc"+Suffix, "RxnCoord"+Suffix});
    PhyQuantName.insert({"mfl"+Suffix, "MeanForceL"+Suffix});
    PhyQuantName.insert({"mfr"+Suffix, "MeanForceR"+Suffix});
    PhyQuantName.insert({"mf"+Suffix, "MeanForce"+Suffix}); 
    OutQuant.insert({"rc"+Suffix, 0.0});
    OutQuant.insert({"mfl"+Suffix, 0.0});
    OutQuant.insert({"mfr"+Suffix, 0.0});
    OutQuant.insert({"mf"+Suffix, 0.0});
}

const prec_t cpihmc::output::internal_energy_estimator() const {return 1.5 * Box.NAtoms * k_B * Box.Temp + Box.QtmKinEng + Box.PotEng;}

void cpihmc::output::update_rxn_coord_quant(const rxn_coord * const RxnCoord, const std::string &Suffix)
{
    OutQuant["rc"+Suffix] = RxnCoord->get_rxn_coord_val();
    OutQuant["mfl"+Suffix] = RxnCoord->get_mean_force_left();
    OutQuant["mfr"+Suffix] = RxnCoord->get_mean_force_right();
    OutQuant["mf"+Suffix] = 0.5 * (OutQuant["mfl"+Suffix]+OutQuant["mfr"+Suffix]);
}

void cpihmc::output::update_physical_quantities(const index_t StepIndex, const vector<rxn_coord *> &RxnCoords)
/* 
    update different physical quantities for each step
    Box: the simulated system
    StepIndex: the index of the current step
*/
{
    OutQuant["ke"] = Box.KinEng;
    OutQuant["pe"] = Box.PotEng;
    OutQuant["etotal"] = Box.TotEng;
    OutQuant["eint"] = internal_energy_estimator();
    if (RxnCoords.size() == 1) update_rxn_coord_quant(RxnCoords[0]);
    else for (index_t i = 0 ; i < RxnCoords.size() ; ++i) update_rxn_coord_quant(RxnCoords[i], "_"+to_string(i));
    OutQuant["ne"] = Box.ElecNum;
    OutQuant["pe_ave"] = (prec_t)StepIndex / (StepIndex+1) * OutQuant["pe_ave"] + Box.PotEng / (StepIndex+1);
}

void cpihmc::output::print_physical_quantities(const index_t StepIndex)
{
    const vector<string> &OutPhyQuant = Input.get_out_phy_quant();
    const unordered_map<string, index_t> &PhyQuantDigits = Input.get_phy_quant_digits();
    OutFile << setw(StepColWidth) << StepIndex;
    for (vector<string>::const_iterator it = OutPhyQuant.begin() ; it < OutPhyQuant.end() ; ++it)
        OutFile << setprecision(PhyQuantDigits.at(*it)) << setw(ColWidth) << OutQuant[*it];
    OutFile << "\n";
}

void cpihmc::output::print_stru_file(const string &FileName, const index_t &NSpin, const bool_t &Direct, const bool_t &Vel,
                                   const bool_t &Magmom) const
{
    const cell &Cell = Box.Cell;
    // ELECTRON_NUMBER
    string Str;
    if (Cell.SetElecNum) Str += "ELECTRON_NUMBER\n" + format("%-.6f\n\n", Cell.ElecNum);
    // ATOMIC_SPECIES
    Str += "ATOMIC_SPECIES\n";
    for(index_t it = 0 ; it < Cell.NType ; ++it)
        Str += format("%s %8.4f %s %s\n", Cell.AtomLabel[it], Cell.AtomMass[it], Cell.PseudoFileName[it], "auto");

    // NUMERICAL_ORBITAL
    if (Cell.OrbitalFileName)
    {
        Str += "\nNUMERICAL_ORBITAL\n";
        for (index_t it = 0 ; it < Cell.NType ; ++it){Str += Cell.OrbitalFileName[it] + "\n";}
    }

    // LATTICE_CONSTANT
    Str += "\nLATTICE_CONSTANT\n" + format("%-.10f\n", Cell.LatConst);

    // LATTICE_VECTORS
    Str += "\nLATTICE_VECTORS\n";
    Str += format("%20.10f%20.10f%20.10f\n", Cell.LatVec.e11, Cell.LatVec.e12, Cell.LatVec.e13);
    Str += format("%20.10f%20.10f%20.10f\n", Cell.LatVec.e21, Cell.LatVec.e22, Cell.LatVec.e23);
    Str += format("%20.10f%20.10f%20.10f\n", Cell.LatVec.e31, Cell.LatVec.e32, Cell.LatVec.e33);

    // ATOMIC_POSITIONS
    Str += "\nATOMIC_POSITIONS\n";
    const string Scale = Direct ? "Direct": "Cartesian";
    Str += Scale + "\n";
    index_t MaskIndex = 0;
    size_t CumAtomNum = 0;
    for (index_t it = 0 ; it < Cell.NType ; ++it)
    {
        Str += "\n" + Cell.Elements[it].Label + " #label\n";
        Str += format("%-8.4f #magnetism\n", Cell.Magnetization[it]);
        Str += format("%d #number of atoms\n", Cell.Elements[it].Number);
        for (index_t ia = 0 ; ia < Cell.Elements[it].Number ; ++ia)
        {
            // output position
            const vec3<prec_t> Coord = Direct ? Cell.Elements[it][ia].Coord / Cell.LatConst * Cell.LatVec.inverse() : Cell.Elements[it][ia].Coord / Cell.LatConst;
            Str += format("%20.10f%20.10f%20.10f", Coord.x, Coord.y, Coord.z);
            Str += format(" m%2d%2d%2d", Cell.Elements[it][ia].Move.x, Cell.Elements[it][ia].Move.y, Cell.Elements[it][ia].Move.z);
            if (Vel) // output velocity
                Str += format(" v%20.10f%20.10f%20.10f", Cell.Elements[it][ia].Vel.x, Cell.Elements[it][ia].Vel.y, Cell.Elements[it][ia].Vel.z);
            if (MaskIndex < Cell.Masks.size())
                if (Cell.Masks[MaskIndex] == CumAtomNum + ia){Str += " mask"; ++MaskIndex;}
            if (NSpin == 2 && Magmom) // output magnetic information
                Str += format(" mag%8.4f", Cell.AtomMagnet[it][ia]);
            else if (NSpin == 4 && Magmom) // output magnetic information
                Str += format(" mag%8.4f%8.4f%8.4f", Cell.AtomMagnetVec[it][ia].x, Cell.AtomMagnetVec[it][ia].y, Cell.AtomMagnetVec[it][ia].z);
            Str += "\n";
        }
        CumAtomNum += Cell.Elements[it].Number;
    }
    ofstream OutFileStru(FileName.c_str());
    OutFileStru << Str;
    OutFileStru.close();
}

void cpihmc::output::print_beads_file(const string &FileName) const
{
    string Str;
    for (index_t i = 0 ; i < Box.NBead ; ++i)
    {
        for (const auto &Bead : Box.Beads) Str += format("%15.10f%15.10f%15.10f\t", Bead.Coords[i].x, Bead.Coords[i].y, Bead.Coords[i].z);
        Str.pop_back();
        Str += "\n";
    }
    ofstream OutFileBeads(FileName.c_str());
    OutFileBeads << Str;
    OutFileBeads.close();
}

void cpihmc::output::save_structure(const index_t StepIndex) const
{
    print_stru_file("ALL_STRU/STRU_"+to_string(StepIndex));
    if (Box.NBead > 1) print_beads_file("ALL_STRU/BEADS_"+to_string(StepIndex));
}

void cpihmc::output::simulation_summary() const
{
    if (Input.get_simu_type() == "HMC" || Input.get_simu_type() == "CHMC")
        cout << "Acceptance Rate:\t" << fixed << setprecision(4) << (Input.get_steps()-RejectTimes) / (prec_t)Input.get_steps() << endl;
    time_t DeltaTime = EndTime - BeginTime;
    index_t Hour = DeltaTime / 3600;
    index_t Minute = (DeltaTime % 3600) / 60;
    index_t Second = DeltaTime % 60;
    printf("Run Time:\t\t%02d:%02d:%02d\n", Hour, Minute, Second);
}