#include "input.h"
using namespace cpihmc;
using namespace std;

bool_t cpihmc::scan_begin(ifstream &InFile, const string &TargetName, const bool_t Restart)
{
    string SearchName;
    bool_t Find = false;
    if (Restart)
    {
        InFile.clear();
        InFile.seekg(0);
    }
    InFile.rdstate();
    while (InFile.good())
    {
        InFile >> SearchName;
        if (SearchName == TargetName)
        {
            Find = true;
            break;
        }
    }
    return Find;
}

cpihmc::input::input(const string &FileName):stru_file("STRU"), pot_type("DP"), deep_pot_model(""), phy_quant_file("PHY_QUANT"), 
                                             col_width(12), simu_type("CHMC"), out_phy_quant({"ke", "pe", "etotal"}), 
                                             phy_quant_digits({{"ke", 6}, {"pe", 6}, {"etotal", 6}, {"eint", 6}, {"pe_ave", 6}, {"pe_var", 6}, 
                                             {"ne", 4}}), steps(10000), phy_quant_intvl(1), stru_intvl(100), temp(300.0), n_type(0), n_evol_step(3), 
                                             time_step(1.0), mass_scal(1.0), hybrid_monte_carlo_ratio(0.8), elec_num_ratio(0.0), mu(0.0), elec_num_range({0.0, 0.0}), 
                                             elec_num_width(0.04), beads_file(""), n_bead(1), n_change_bead(1), model_devi_intvl(100), model_devi_file("MODEL_DEVI")
{read(FileName);}

void cpihmc::input::insert_rxn_coord_digits(const std::string &Suffix)
{
    phy_quant_digits.insert({"rc"+Suffix, 2});
    phy_quant_digits.insert({"mfl"+Suffix, 6});
    phy_quant_digits.insert({"mfr"+Suffix, 6});
    phy_quant_digits.insert({"mf"+Suffix, 6});
}

void cpihmc::input::read(const string &FileName)
{
    ifstream InFileParam(FileName.c_str(), ios::in);
    if (!InFileParam) throw std::invalid_argument("Can't find the INPUT file.");
    InFileParam.clear();
    InFileParam.seekg(0); // back to position 0

    char Word[80];
    bool_t FileErr = true;
    while (InFileParam.good())
    {
        InFileParam >> Word;
        InFileParam.ignore(150, '\n');
        if (strcmp(Word , "INPUT_PARAMETERS") == 0)
        {
            FileErr = false;
            break;
        }
    }
    if (FileErr) throw std::invalid_argument("Error parameter list.");

    char WordLower[80];
    while (InFileParam.good())
    {
        InFileParam >> Word;
        str2lower(Word, WordLower);
        if (strcmp("stru_file", WordLower) == 0)
            read_value(InFileParam, stru_file);
        else if (strcmp("pot_type", WordLower) == 0)
            read_value(InFileParam, pot_type);
        else if (strcmp("deep_pot_model", WordLower) == 0)
            read_value(InFileParam, deep_pot_model);
        else if (strcmp("simu_type", WordLower) == 0)
            read_value(InFileParam, simu_type);
        else if (strcmp("steps", WordLower) == 0)
            read_value(InFileParam, steps);
        else if (strcmp("temp", WordLower) == 0)
            read_value(InFileParam, temp);
        else if (strcmp("n_type", WordLower) == 0)
            read_value(InFileParam, n_type);
        else if (strcmp("element_type", WordLower) == 0)
        {
            element_type.resize(n_type);
            for (index_t i = 0 ; i < n_type - 1 ; ++i)
                InFileParam >> element_type[i];
            read_value(InFileParam, element_type[n_type-1]);
        }
        else if (strcmp("element_index", WordLower) == 0)
        {
            read_value(InFileParam, element_index);
            assert(n_type == element_index.size());
        }
        else if (strcmp("wall", WordLower) == 0)
        {
            string Name;
            prec_t Pos;
            while (InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#')
                {
                    InFileParam.ignore(150, '\n');
                    break;
                }
                else
                {
                    InFileParam >> Name >> Pos;
                    assert(WallIndexes.find(Name) != WallIndexes.end());
                    assert(wall.find(Name) == wall.end());
                    wall.insert({Name, Pos});
                }
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
        }
        else if (strcmp("group", WordLower) == 0)
        {
            string Name;
            array<index_t, 3> Range;
            vector<index_t> Indexes;
            while (InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#'){InFileParam.ignore(150, '\n'); break;}
                else
                {
                    InFileParam >> Name;
                    assert(group_indexes.find(Name) == group_indexes.end());
                    Range[2] = 1;
                    char NextChar;
                    string Num;
                    index_t Index = 0;
                    do
                    {
                        if (NextChar != ':') Num += InFileParam.get();
                        else
                        {
                            InFileParam.get();
                            Range[Index] = stoi(Num);
                            ++Index;
                            Num = "";
                        }
                        NextChar = InFileParam.peek();
                    } while (NextChar != '\n' && NextChar != ' ' && NextChar != '\t' && NextChar != '#');
                    Range[Index] = stoi(Num);
                    if (Index != 1 && Index != 2) throw invalid_argument("Only two or three integers are allowed to defined the group.");
                    assert(Range[2]);
                    if ((Range[1]-Range[0])*Range[2]<0) throw invalid_argument("The group range of the parameter `Group` is invalid.");
                    const size_t IndexNum = (Range[1]-Range[0])/Range[2];
                    Indexes.resize(IndexNum);
                    for (index_t i = 0 ; i < IndexNum ; ++i) Indexes[i] = Range[0] + i * Range[2];
                    group_indexes.insert({Name, Indexes});
                }
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
        }
        else if (strcmp("n_evol_step", WordLower) == 0)
            read_value(InFileParam, n_evol_step);
        else if (strcmp("time_step", WordLower) == 0)
            read_value(InFileParam, time_step);
	    else if (strcmp("mass_scal", WordLower) == 0)
	        read_value(InFileParam, mass_scal);
        else if (strcmp("hybrid_monte_carlo_ratio", WordLower) == 0)
            read_value(InFileParam, hybrid_monte_carlo_ratio);
        else if (strcmp("virt_atom", WordLower) == 0)
        {
            string Name;
            index_t VirtAtomIndex = 0;
            index_t Index;
            vector<index_t> Indexes;
            while (InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#')
                {
                    InFileParam.ignore(150, '\n');
                    break;
                }
                else if ((InFileParam.peek() >= 'a' && InFileParam.peek() <= 'z') || (InFileParam.peek() >= 'A' && InFileParam.peek() <= 'Z'))
                {
                    if (Name.size())
                    {
                        assert(Indexes.size() > 1);
                        virt_atom.insert(pair<string, index_t>{Name, VirtAtomIndex});
                        virt_atom_names.push_back(Name);
                        ++VirtAtomIndex;
                        virt_atom_indexes.push_back(Indexes);
                    }
                    InFileParam >> Name;
                    Indexes.clear();
                }
                else
                {
                    InFileParam >> Index;
                    assert(Index >= 0);
                    Indexes.push_back(Index);
                }
            }
            if (Name.size())
            {
                assert(Indexes.size() > 1);
                virt_atom.insert(pair<string, index_t>{Name, VirtAtomIndex});
                virt_atom_names.push_back(Name);
                virt_atom_indexes.push_back(Indexes);
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
        }
        else if (strcmp("rxn_coord", WordLower) == 0)
        {
            string Type;
            string Index;
            vector<index_t> AtomIndexes;
            string Val;
            prec_t Param;
            vector<prec_t> Params;
            while (InFileParam.good() && InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#')
                {
                    InFileParam.ignore(150, '\n');
                    break;
                }
                else
                {
                    InFileParam >> Type;
                    assert(RxnCoordAtomNums.find(Type) != RxnCoordAtomNums.end());
                    AtomIndexes.clear();
                    Params.clear();
                    for (index_t i = 0 ; i < RxnCoordAtomNums.at(Type) ; ++i)
                    {
                        InFileParam >> Index;
                        if ((Index[0] >= 'a' && Index[0] <= 'z') || (Index[0] >= 'A' && Index[0] <= 'Z')) AtomIndexes.push_back(-1-virt_atom[Index]);
                        else AtomIndexes.push_back(stoi(Index));
                    }
                    InFileParam >> Val;
                    for (index_t i = 0 ; i < RxnCoordParamNums.at(Type) ; ++i)
                    {
                        InFileParam >> Param;
                        Params.push_back(Param);
                    }
                    if (strcmp(Val.c_str(), "FILE") == 0) rxn_coord.push_back(rxn_coord_info{Type, AtomIndexes, 0.0, Params, false});
                    else rxn_coord.push_back(rxn_coord_info{Type, AtomIndexes, (prec_t)stod(Val), Params});
                }
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
            if (rxn_coord.size() == 1) insert_rxn_coord_digits();
            else for (index_t i = 0 ; i < rxn_coord.size() ; ++i) insert_rxn_coord_digits("_"+to_string(i));
        }
        else if (strcmp("elec_num_ratio", WordLower) == 0)
            read_value(InFileParam, elec_num_ratio);
        else if (strcmp("mu", WordLower) == 0)
            read_value(InFileParam, mu);
        else if (strcmp("elec_num_range", WordLower) == 0)
        {
            InFileParam >> elec_num_range.first >> elec_num_range.second;
            string Line;
            getline(InFileParam, Line);
        }
        else if (strcmp("elec_num_width", WordLower) == 0)
            read_value(InFileParam, elec_num_width);
        else if (strcmp("beads_file", WordLower) == 0)
            read_value(InFileParam, beads_file);
        else if (strcmp("n_bead", WordLower) == 0)
            read_value(InFileParam, n_bead);
        else if (strcmp("n_change_bead", WordLower) == 0)
            read_value(InFileParam, n_change_bead);
        else if (strcmp("bead_index", WordLower) == 0)
        {
            string Name;
            index_t Index;
            while (InFileParam.good() && InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#'){InFileParam.ignore(150, '\n'); break;}
                else if ((InFileParam.peek() >= 'a' && InFileParam.peek() <= 'z') || (InFileParam.peek() >= 'A' && InFileParam.peek() <= 'Z'))
                {
                    InFileParam >> Name;
                    bead_index.insert(bead_index.end(), group_indexes[Name].begin(), group_indexes[Name].end());
                }
                else
                {
                    InFileParam >> Index;
                    bead_index.push_back(Index);
                }
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
            for (const auto &Index : bead_index) assert(Index >= 0);
            for (index_t i = 0 ; i < bead_index.size() ; ++i)
                for (index_t j = i + 1 ; j < bead_index.size() ; ++j)
                    if (bead_index[i] == bead_index[j]) throw invalid_argument("Duplicate indexes aren't allowed.");
        }
        else if (strcmp("model_devi_deep_pot_models", WordLower) == 0)
            read_value(InFileParam, model_devi_deep_pot_models);
        else if (strcmp("model_devi_intvl", WordLower) == 0)
            read_value(InFileParam, model_devi_intvl);
        else if (strcmp("model_devi_file", WordLower) == 0)
            read_value(InFileParam, model_devi_file);
        else if (strcmp("phy_quant_file", WordLower) == 0)
            read_value(InFileParam, phy_quant_file);
        else if (strcmp("col_width", WordLower) == 0)
            read_value(InFileParam, col_width);
        else if (strcmp("out_phy_quant", WordLower) == 0)
        {
            out_phy_quant.clear();
            read_value(InFileParam, out_phy_quant);
        }
        else if (strcmp("phy_quant_digits", WordLower) == 0)
        {
            string PhyQuant;
            index_t Digit;
            while (InFileParam.good() && InFileParam.peek() != '\n')
            {
                if (InFileParam.peek() == ' ' || InFileParam.peek() == '\t') InFileParam.get();
                else if (InFileParam.peek() == '#')
                {
                    InFileParam.ignore(150, '\n');
                    break;
                }
                else
                {
                    InFileParam >> PhyQuant >> Digit;
                    assert(Digit != 0);
                    assert(phy_quant_digits.find(PhyQuant) != phy_quant_digits.end());
                    phy_quant_digits[PhyQuant] = Digit;
                }
            }
            if (InFileParam.peek() == '\n') InFileParam.get(); // get rid of '\n'
        }
        else if (strcmp("phy_quant_intvl", WordLower) == 0)
            read_value(InFileParam, phy_quant_intvl);
        else if (strcmp("stru_intvl", WordLower) == 0)
            read_value(InFileParam, stru_intvl);
        else
        {
            if (WordLower[0] != '#' && WordLower[0] != '/')
                cerr << " THE PARAMETER NAME '" << Word << "' IS NOT USED!" << endl;
            InFileParam.ignore(150, '\n');
        }
        if (InFileParam.peek() == EOF || InFileParam.eof()) break;
        else if (InFileParam.bad()) cerr << " Bad input parameters. " << endl;
        else if (InFileParam.fail())
        {
            cerr << " word = " << Word << "\n" << " Fail to read parameters. " << endl;
            InFileParam.clear();
        }
        else if (!InFileParam.good()) break;
    }
}

void cpihmc::input::print(const string &FileName) const
{
    ofstream OutFileParam(FileName.c_str());
    OutFileParam << "INPUT_PARAMETERS" << endl;
    OutFileParam << setiosflags(ios::left);
    OutFileParam << "#Parameters (General)" << endl;
    out_param(OutFileParam, "Stru_File", stru_file, "The structure file.");
    out_param(OutFileParam, "Pot_Type", pot_type, "The type of potential energy and force calculator.");
    out_param(OutFileParam, "Deep_Pot_Model", deep_pot_model, "The DP potential for the system.");
    out_param(OutFileParam, "Simu_Type", simu_type, "The type of the simulation.");
    out_param(OutFileParam, "Steps", steps, "The number of the simulation steps.");
    OutFileParam << "#Parameters (System)" << endl;
    out_param(OutFileParam, "Temp", temp, "The temperature of the system.");
    out_param(OutFileParam, "N_Type", n_type, "The number of element(s).");
    stringstream StrSET;
    for (index_t i = 0 ; i < n_type ; ++i) StrSET << element_type[i] << " ";
    string StrET = StrSET.str();
    StrET.pop_back();
    out_param(OutFileParam, "Element_Type", StrET, "Different element type(s) at the order in DP potential when the paramter \"Element_Index\" is not set.");
    if (element_index.size())
    {
        stringstream StrSEI;
        for (index_t i = 0 ; i < n_type ; ++i) StrSEI << element_index[i] << " ";
        string StrEI = StrSEI.str();
        StrEI.pop_back();
        out_param(OutFileParam, "Element_Index", StrEI, "Different element index(s) used in DP potential mapping to element type(s) one by one.");
    }
    if (wall.size())
    {
        stringstream StrSW;
        for (const auto &Wall : wall) StrSW << Wall.first << " " << Wall.second << " ";
        string StrW = StrSW.str();
        StrW.pop_back();
        out_param(OutFileParam, "Wall", StrW, "The wall(s) set to limit atom coordinates in a cell.");
    }
    if (simu_type == "MD")
    {
        OutFileParam << "#Parameters (Molecular Dynamics)" << endl;
        out_param(OutFileParam, "Time_Step", time_step, "The time step of Molecular Dynamics simulation.");
    }
    else if (simu_type == "HMC" || simu_type == "CHMC")
    {
        OutFileParam << "#Parameters (Hybrid Monte Carlo)" << endl;
        out_param(OutFileParam, "N_Evol_Step", n_evol_step, "The number of evolution steps between two judges in the HMC scheme.");
        out_param(OutFileParam, "Time_Step", time_step, "The time step of the HMC scheme.");
        out_param(OutFileParam, "Mass_Scal", mass_scal, "The mass scaling coefficient used in the HMC scheme.");
        if (simu_type == "CHMC")
        {
            OutFileParam << "#Parameters (Constraint)" << endl;
            out_param(OutFileParam, "Hybrid_Monte_Carlo_Ratio", hybrid_monte_carlo_ratio, "The ratio of the HMC part in the centroid move.");
            if (virt_atom.size())
            {
                stringstream StrSVA;
                for (unordered_map<string, index_t>::const_iterator IterVirtAtom = virt_atom.begin() ; IterVirtAtom != virt_atom.end() ; ++IterVirtAtom)
                {
                    StrSVA << IterVirtAtom->first << " ";
                    for (const auto &Index : virt_atom_indexes[IterVirtAtom->second]) StrSVA << Index << " ";
                }
                string StrVA = StrSVA.str();
                StrVA.pop_back();
                out_param(OutFileParam, "Virt_Atom", StrVA, "The virtual atom(s) used in the constraint.");
            }
            stringstream StrSRC;
            for (const auto &RxnCoord : rxn_coord)
            {
                StrSRC << RxnCoord.type << " ";
                for (const auto &AtomIndex : RxnCoord.atom_indexes)
                {
                    if (AtomIndex < 0) StrSRC << virt_atom_names[-1-AtomIndex] << " ";
                    else StrSRC << AtomIndex << " ";
                }
                if (RxnCoord.set_rxn_coord) StrSRC << RxnCoord.value << " ";
                else StrSRC << "FILE ";
                for (const auto &Param : RxnCoord.params) StrSRC << Param << " ";
            }
            string StrRC = StrSRC.str();
            StrRC.pop_back();
            out_param(OutFileParam, "Rxn_Coord", StrRC, "The reaction coordinate(s) used as the constraint.");
            if (elec_num_ratio)
            {
                OutFileParam << "#Parameters (Grand Canonical Ensemble)" << endl;
                out_param(OutFileParam, "Elec_Num_Ratio", elec_num_ratio, "The ratio of the electron number change part in all moves.");
                out_param(OutFileParam, "Mu", mu, "The parameter `mu` for the grand canonical ensemble.");
                stringstream StrSENR;
                StrSENR << elec_num_range.first << " " << elec_num_range.second;
                out_param(OutFileParam, "Elec_Num_Range", StrSENR.str(), "The range of the electron number.");
                out_param(OutFileParam, "Elec_Num_Width", elec_num_width, "The change width of the electron number.");
            }
        }
    }
    if (n_bead > 1)
    {
        OutFileParam << "#Parameters (Path Integral)" << endl;
        out_param(OutFileParam, "Beads_File", beads_file, "The initial beads coordinates file.");
        out_param(OutFileParam, "N_Bead", n_bead, "The number of beads.");
        out_param(OutFileParam, "N_Change_Bead", n_change_bead, "The number of bead(s) changed in each step.");
        stringstream StrSBI;
        for (index_t i = 0 ; i < bead_index.size() ; ++i) StrSBI << bead_index[i] << " ";
        string StrBI = StrSBI.str();
        StrBI.pop_back();
        out_param(OutFileParam, "Bead_Index", StrBI, "The index(s) of different bead(s) in coordinates.");
    }
    if (model_devi_deep_pot_models.size())
    {
        OutFileParam << "#Parameters (Model Deviation)" << endl;
        stringstream StrSMDDPM;
        for (const auto &DeepPotModel : model_devi_deep_pot_models) StrSMDDPM << DeepPotModel << " ";
        string StrMDDPM = StrSMDDPM.str();
        StrMDDPM.pop_back();
        out_param(OutFileParam, "Model_Devi_Deep_Pot_Models", StrMDDPM, "The DP models used for the model deviation.");
        out_param(OutFileParam, "Model_Devi_Intvl", model_devi_intvl, "The interval to output the model deviation results.");
        out_param(OutFileParam, "Model_Devi_File", model_devi_file, "The file to output the model deviation results.");
    }
    OutFileParam << "#Parameters (Output)" << endl;
    out_param(OutFileParam, "Phy_Quant_File", phy_quant_file, "The file to output the result.");
    out_param(OutFileParam, "Col_Width", col_width, "The width of the column in the file to output the result.");
    stringstream StrSOPQ;
    for (const auto &PhyQuant : out_phy_quant) StrSOPQ << PhyQuant << " ";
    string StrOPQ = StrSOPQ.str();
    StrOPQ.pop_back();
    out_param(OutFileParam, "Out_Phy_Quant", StrOPQ, "The physical quantity(ies) to be output.");
    stringstream StrSPQD;
    for (const auto &PhyQuant : out_phy_quant) StrSPQD << PhyQuant << " " << phy_quant_digits.at(PhyQuant) << " ";
    string StrPQD = StrSPQD.str();
    StrPQD.pop_back();
    out_param(OutFileParam, "Phy_Quant_Digits", StrPQD, "The decimal digits of physical quantity(ies) to be output.");
    out_param(OutFileParam, "Phy_Quant_Intvl", phy_quant_intvl, "The number of steps in the interval to output the physical quantity(ies).");
    out_param(OutFileParam, "Stru_Intvl", stru_intvl, "The number of steps in the interval to output the structure.");
    OutFileParam.close();
}