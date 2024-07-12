//
// Created by Jin Bin on 2021/11/03.
//

#include "constant.h"

constant::constant():ntype(0), T(300), M_scaling(1), N_MAX(1293), dump_steps(100), steps(100000), Delta_t(2.4), reaction_coordinate(0), 
hmc_proportion(0.5), N_step(3), radial(0.4), length_step(0.3), theta_step(0.1309), phi_step(1.0472), work_path("."),
atom_file("STRU"), out_file("energy.dat"), initial_beads(""), beads_file("BEADS"), potential_type("DP"), DP_model("graph.pb"),
calculation_type("CHMC"), RC_type("Difference"), P(1), J(0), bead_num(1), output_terms({"ke", "pe", "etotal", "rc", "mfl", "mfr"}),
quantity_digits({{"ke", 6}, {"pe", 6}, {"etotal", 6}, {"eint", 6}, {"rc", 2}, {"mfl", 6}, {"mfr", 6}, {"mf", 6}, {"pe_ave", 6}, {"pe_var", 6}, {"ne", 2}}),
DP_model_prefix(""), Ne_range(pair<double, double>{0.5, 0.66}), mu(-0.5), delta_Ne(0.04), GC_type("Continuous"), Ne_proportion(0.1),
model_devi_interval(100), model_devi_file("model_devi.out"), set_rc(false){
    bead_index.resize(bead_num);
    for (int i = 0 ; i < bead_num ; ++i){bead_index[i]=0;}
    constraint_atom = new int [3];
    for (int i = 0 ; i < 3 ; ++i){constraint_atom[i]=0;}
}

constant::constant(const string &input):ntype(0), T(300), M_scaling(1), N_MAX(1293), dump_steps(100), steps(100000), Delta_t(2.4), reaction_coordinate(0), 
hmc_proportion(0.5), N_step(3), radial(0.4), length_step(0.3), theta_step(0.1309), phi_step(1.0472), work_path("."),
atom_file("STRU"), out_file("energy.dat"), initial_beads(""), beads_file("BEADS"), potential_type("DP"), DP_model("graph.pb"), 
calculation_type("CHMC"), RC_type("Difference"), P(1), J(0), bead_num(1), output_terms({"ke", "pe", "etotal", "rc", "mfl", "mfr"}),
quantity_digits({{"ke", 6}, {"pe", 6}, {"etotal", 6}, {"eint", 6}, {"rc", 2}, {"mfl", 6}, {"mfr", 6}, {"mf", 6}, {"pe_ave", 6}, {"pe_var", 6}, {"ne", 2}}),
DP_model_prefix(""), Ne_range(pair<double, double>{0.5, 0.66}), mu(-0.5), delta_Ne(0.04), GC_type("Continuous"), Ne_proportion(0.1),
model_devi_interval(100), model_devi_file("model_devi.out"), set_rc(false){
    bead_index.resize(bead_num);
    for (int i = 0 ; i < bead_num ; ++i){bead_index[i]=0;}
    constraint_atom = new int [3];
    for (int i = 0 ; i < 3 ; ++i){constraint_atom[i]=0;}
    Read(input);
}

constant::~constant(){
    delete [] constraint_atom;
}

void constant::Read(const string &input){
    ifstream ifs(input.c_str(), ios::in);
    if (!ifs)
    {
        cout << " Can't find the INPUT file." << endl;
        return;
    }
    ifs.clear();
    ifs.seekg(0); // back to position 0

    char word[80];
    char word1[80];
    int ierr = 0;

    while (ifs.good()){
        ifs >> word;
        ifs.ignore(150, '\n');
        if (strcmp(word , "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
    }

    if (ierr == 0)
        cout << " Error parameter list." << endl;

    while (ifs.good()){
        ifs >> word1;
        if (ifs.eof()) break;
        strtolower(word1, word);
        if (strcmp("ntype", word) == 0)
            read_value(ifs, ntype);
        else if (strcmp("t", word) == 0)
            read_value(ifs, T);
	    else if (strcmp("m_scaling", word) == 0)
	        read_value(ifs, M_scaling);
        else if (strcmp("n_max", word) == 0)
            read_value(ifs, N_MAX);
        else if (strcmp("dump", word) == 0)
            read_value(ifs, dump_steps);
        else if (strcmp("interval", word) == 0)
            read_value(ifs, model_devi_interval);
        else if (strcmp("steps", word) == 0)
            read_value(ifs, steps);
        else if (strcmp("output", word) == 0){
            output_terms.clear();
            read_value(ifs, output_terms);
        }
        else if (strcmp("digits", word) == 0){
            string phys_quan;
            int digit;
            while (ifs.peek() != '\n')
            {
                if (ifs.peek() == ' ' || ifs.peek() == '\t')
                    ifs.get();
                else if (ifs.peek() == '#'){
                    ifs.ignore(150, '\n');
                    break;
                }
                else {
                    ifs >> phys_quan >> digit;
                    assert(digit != 0);
                    assert(quantity_digits.find(phys_quan)!=quantity_digits.end());
                    quantity_digits[phys_quan] = digit;
                }
            }
            if (ifs.peek() == '\n')
                ifs.get(); // get rid of '\n'
        }
        else if (strcmp("model_devi_models", word) == 0)
            read_value(ifs, model_devi_DP_models);
        else if (strcmp("delta_t", word) == 0)
            read_value(ifs, Delta_t);
        else if (strcmp("element_type", word) == 0){
            element_type.resize(ntype);
            string temp;
            for (int i = 0 ; i < ntype - 1 ; ++i){
                ifs >> temp;
                element_type[i] = temp;
            }
            read_value(ifs, temp);
            element_type[ntype-1] = temp;
        }
        else if (strcmp("element_index", word) == 0){
            read_value(ifs, element_index);
            assert(ntype == element_index.size());
        }
        else if (strcmp("wall", word) == 0)
            read_value(ifs, wall);
        else if (strcmp("w_pos", word) == 0){
            read_value(ifs, wall_pos);
            assert(wall.size() == wall_pos.size());
        }
        else if (strcmp("rc_type", word) == 0)
            read_value(ifs, RC_type);
        else if (strcmp("constraint_atom", word) == 0){
            if (strcmp(RC_type.c_str(), "Difference") == 0){
                ifs >> constraint_atom[0] >> constraint_atom[1];
                read_value(ifs, constraint_atom[2]);
            }
            else if (strcmp(RC_type.c_str(), "Distance") == 0){
                delete [] constraint_atom;
                constraint_atom = new int [2];
                ifs >> constraint_atom[0];
                read_value(ifs, constraint_atom[1]);
            }
        }
        else if (strcmp("calculation", word) == 0)
            read_value(ifs, calculation_type);
        else if (strcmp("potential", word) == 0)
            read_value(ifs, potential_type);
        else if (strcmp("reaction_coordinate", word) == 0){
            read_value(ifs, reaction_coordinate);
            set_rc = true;
        }
        else if (strcmp("hmc_proportion", word) == 0)
            read_value(ifs, hmc_proportion);
        else if (strcmp("n_step", word) == 0)
            read_value(ifs, N_step);
        else if (strcmp("mu", word) == 0)
            read_value(ifs, mu);
        else if (strcmp("delta_ne", word) == 0)
            read_value(ifs, delta_Ne);
        else if (strcmp("gc_type", word) == 0)
            read_value(ifs, GC_type);
        else if (strcmp("ne_proportion", word) == 0)
            read_value(ifs, Ne_proportion);
        else if (strcmp("ne_range", word) == 0){
            double low, high;
            ifs >> low;
            read_value(ifs, high);
            Ne_range = pair<double, double>{low, high};
        }
        else if (strcmp("radial", word) == 0)
            read_value(ifs, radial);
        else if (strcmp("hmc_lstep", word) == 0)
            read_value(ifs, length_step);
        else if (strcmp("hmc_theta_step", word) == 0)
            read_value(ifs, theta_step);
        else if (strcmp("hmc_phi_step", word) == 0)
            read_value(ifs, phi_step);
        else if (strcmp("mc_theta_step", word) == 0)
            read_value(ifs, theta_step);
        else if (strcmp("mc_phi_step", word) == 0)
            read_value(ifs, phi_step);
        else if (strcmp("work_path", word) == 0)
            read_value(ifs, work_path);
        else if (strcmp("p", word) == 0)
            read_value(ifs, P);
        else if (strcmp("j", word) == 0)
            read_value(ifs, J);
        else if (strcmp("n_bead", word) == 0)
            read_value(ifs, bead_num);
        else if (strcmp("bead_index", word) == 0){
            bead_index.resize(bead_num);
            int temp;
            for (int i = 0 ; i < bead_num - 1 ; ++i){
                ifs >> temp;
                bead_index[i] = temp;
            }
            read_value(ifs, temp);
            bead_index[bead_num-1] = temp;
        }
        else if (strcmp("atom_file", word) == 0)
            read_value(ifs, atom_file);
        else if (strcmp("out_file", word) == 0)
            read_value(ifs, out_file);
        else if (strcmp("model_devi_file", word) == 0)
            read_value(ifs, model_devi_file);
        else if (strcmp("initial_beads", word) == 0)
            read_value(ifs, initial_beads);
        else if (strcmp("beads_file", word) == 0)
            read_value(ifs, beads_file);
        else if (strcmp("dp_model", word) == 0)
            read_value(ifs, DP_model);
        else if (strcmp("dp_model_prefix", word) == 0)
            read_value(ifs, DP_model_prefix);
        else {
            if(word[0] != '#' && word[0] != '/')
                cout<<" THE PARAMETER NAME '" << word << "' IS NOT USED!" << endl;
            ifs.ignore(150, '\n');
        }
        if (ifs.eof())
            break;
        else if (ifs.bad())
            cout << " Bad input parameters. " << endl;
        else if (ifs.fail()) {
            cout << " word = " << word << endl;
            cout << " Fail to read parameters. " << endl;
            ifs.clear();
        }
        else if (ifs.good() == 0)
            break;
    }
}

void constant::Print(const string &output) const{
    char cal_type[16];
    constant::strtolower(calculation_type.c_str(), cal_type);

    ofstream ofs(output.c_str());
    ofs << "INPUT_PARAMETERS" << endl;
    ofs << setiosflags(ios::left);
    ofs << "#Parameters (General)" << endl;
    OUTP(ofs, "DP_model", DP_model, "The DP potential for the system.");
    OUTP(ofs, "work_path", work_path, "The work path for ab initio calculation.");
    OUTP(ofs, "atom_file", atom_file, "The structure file.");
    OUTP(ofs, "out_file", out_file, "The file to output the result.");
    OUTP(ofs, "calculation", calculation_type, "The type of calculation.");
    OUTP(ofs, "potential", potential_type, "The type of potential energy and force calculator.");
    stringstream so;
    int n_output = output_terms.size();
    for (int i = 0 ; i < n_output-1 ; ++i)
        so << output_terms[i] << " ";
    so << output_terms[n_output-1];
    OUTP(ofs, "OUTPUT", so.str(), "The physical quantity(ies) to be output.");
    stringstream sq;
    for (int i = 0 ; i < n_output-1 ; ++i)
        sq << output_terms[i] << " " << quantity_digits.at(output_terms[i]) << " ";
    sq << output_terms[n_output-1] << " " << quantity_digits.at(output_terms[n_output-1]);
    OUTP(ofs, "digits", sq.str(), "The decimal digits of physical quantity(ies) to be output.");
    OUTP(ofs, "steps", steps, "The number of the simulation steps.");
    OUTP(ofs, "dump", dump_steps, "The number of steps in the interval to output the structure.");
    ofs << "#Parameters (System)" << endl;
    OUTP(ofs, "T", T, "The temperature of the simulation.");
    OUTP(ofs, "ntype", ntype, "The number of element(s).");
    stringstream se;
    for (int i = 0 ; i < ntype - 1 ; ++i)
        se << element_type[i] << " ";
    se << element_type[ntype-1];
    OUTP(ofs, "element_type", se.str(), "Different element type(s) at the order in DP potential when the paramter \"element_index\" is not set.");
    if (element_index.size()){
        stringstream sei;
        for (int i = 0 ; i < ntype - 1 ; ++i)
            sei << element_index[i] << " ";
        sei << element_index[ntype-1];
        OUTP(ofs, "element_index", sei.str(), "Different element index(s) used in DP potential mapping to element type(s) one by one.");
    }
    if (strcmp(cal_type, "chmc") == 0 || strcmp(cal_type, "cmc") == 0){
        OUTP(ofs, "RC_type", RC_type, "The type of reaction coordinate");
        OUTP(ofs, "reaction_coordinate", reaction_coordinate, "The reaction coordinate value.");
        stringstream ss;
        if (strcmp(RC_type.c_str(), "Difference") == 0)
            ss << constraint_atom[0] << " " << constraint_atom[1] << " " << constraint_atom[2];
        else if (strcmp(RC_type.c_str(), "Distance") == 0)
            ss << constraint_atom[0] << " " << constraint_atom[1];
        OUTP(ofs, "constraint_atom", ss.str(), "The atom index about the reaction coordinate.");
    }
    int n_wall = wall.size();
    if (n_wall) {
        stringstream sw;
        for (int i = 0 ; i < n_wall-1 ; ++i)
            sw << wall[i] << " ";
        sw << wall[n_wall-1];
        OUTP(ofs, "Wall", sw.str(), "The wall(s) set to limit atom coordinates in a cell.");
        stringstream swp;
        vector<double>::const_iterator it = wall_pos.begin();
        for ( ; it < wall_pos.end() - 1 ; ++it)
            swp << *it << " ";
        swp << *it;
        OUTP(ofs, "W_pos", swp.str(), "The position(s) of the set wall(s).");
    }
    else
        OUTP(ofs, "Wall", "None", "The wall(s) set to limit atom coordinates in a cell.");
    OUTP(ofs, "N_MAX", N_MAX, "Maximum number of atoms in the structure.");
    if (strlen(DP_model_prefix.c_str())){
        ofs << "#Parameters (Grand Canonical Ensemble)" << endl;
        OUTP(ofs, "DP_model_prefix", DP_model_prefix, "The prefix of the DP potential with different electron numbers or the DP potential containing electron number.");
        OUTP(ofs, "GC_type", GC_type, "The type of the DP potential. Options: Continuous, Discrete.");
        stringstream sne;
        sne << Ne_range.first << " " << Ne_range.second;
        OUTP(ofs, "Ne_range", sne.str(), "The range of the electron number.");
        OUTP(ofs, "delta_Ne", delta_Ne, "The change interval of the electron number.");
        OUTP(ofs, "Ne_proportion", Ne_proportion, "The proportion of the electron number change part in all moves.");
        OUTP(ofs, "mu", mu, "The parameter mu for the grand canonical ensemble.");
    }
    if (model_devi_DP_models.size()){
        ofs << "#Parameters (Model Deviation)" << endl;
        stringstream smdm;
        vector<string>::const_iterator it = model_devi_DP_models.begin();
        for ( ; it < model_devi_DP_models.end() - 1 ; ++it)
            smdm << *it << " ";
        smdm << *it;
        OUTP(ofs, "model_devi_models", smdm.str(), "The DP models used for the model deviation.");
        OUTP(ofs, "model_devi_file", model_devi_file, "The file to output the model deviation results.");
        OUTP(ofs, "interval", model_devi_interval, "The interval to output the model deviation results.");
    }
    if (strcmp(cal_type, "md") == 0){
        ofs << "#Parameters (Molecular Dynamics)" << endl;
        OUTP(ofs, "Delta_t", Delta_t, "The time interval of Molecular Dynamics simulation.");
    }
    if (strcmp(cal_type, "hmc") == 0){
        ofs << "#Parameters (Hamiltonian Monte Carlo)" << endl;
        OUTP(ofs, "M_scaling", M_scaling, "The mass scaling coefficient used in the HMC scheme.");
        OUTP(ofs, "Delta_t", Delta_t, "The time interval of Hamiltonian Monte Carlo simulation.");
        OUTP(ofs, "N_step", N_step, "The number of steps between two judge.");
    }
    if (strcmp(cal_type, "chmc") == 0){
        ofs << "#Parameters (Constrained Hamiltonian Monte Carlo)" << endl;
        OUTP(ofs, "M_scaling", M_scaling, "The mass scaling coefficient used in the HMC scheme.");
        OUTP(ofs, "Delta_t", Delta_t, "The time interval of Hamiltonian Monte Carlo simulation.");
        OUTP(ofs, "HMC_proportion", hmc_proportion, "The proportion of the HMC part in centroid move.");
        OUTP(ofs, "N_step", N_step, "The number of steps between two judge.");
        if (strcmp(RC_type.c_str(), "Difference") == 0){
            OUTP(ofs, "HMC_lstep", length_step, "The length step of MC movement in HMC method.");
            OUTP(ofs, "HMC_theta_step", theta_step, "The theta step of MC movement in HMC method.");
            OUTP(ofs, "HMC_phi_step", phi_step, "The phi step of MC movement in HMC method.");
        }
        else if (strcmp(RC_type.c_str(), "Distance") == 0)
        {
            OUTP(ofs, "radial", radial, "The radial of ball to adjust angle.");
        }
    }
    if (strcmp(cal_type, "cmc") == 0){
        ofs << "#Parameters (Constrained Monte Carlo)" << endl;
        OUTP(ofs, "radial", radial, "The radial of ball to adjust angle.");
        if (strcmp(RC_type.c_str(), "Difference") == 0){
            OUTP(ofs, "MC_theta_step", theta_step, "The theta step of central atom in the constraint.");
            OUTP(ofs, "MC_phi_step", phi_step, "The phi step of central atom in the constraint.");
        }
    }
    if (P > 1){
        ofs << "#Parameters (Path Integral)" << endl;
        OUTP(ofs, "initial_beads", initial_beads, "The initial beads coordinates file.");
        OUTP(ofs, "beads_file", beads_file, "The prefix of beads coordinates file.");
        OUTP(ofs, "P", P, "The number of beads.");
        OUTP(ofs, "J", J, "The number of beads changed in each step.");
        OUTP(ofs, "N_bead", bead_num, "The number of bead atom(s).");
        stringstream sb;
        for (int i = 0 ; i < bead_num - 1 ; ++i)
            sb << bead_index[i] << " ";
        sb << bead_index[bead_num-1];
        OUTP(ofs, "bead_index", sb.str(), "The index(s) of different bead(s) in coordinates.");
    }
    ofs.close();
}

void constant::strtolower(const char *sa, char *sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; ++i)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}
