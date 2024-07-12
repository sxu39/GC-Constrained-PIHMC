#include <iostream>
#include "cell.h"
#include "MD_NVE.h"
#include "HMC_muVT.h"
#include "run_time.h"
#include "ml_potential.h"
#include "first_principle.h"
#include "gen_func.h"
#include "model_devi.h"

// define these variables as global variables to support its usage in other source files via extern
constant Consts("INPUT"); // obtain parameters
string atom_file = Consts.work_path + "/" + Consts.atom_file; // structure file
cell Cell(atom_file); // obtain cell from structure file, class cell is based on corresponding part of ABACUS

int digits = floor(log10(max(1, Consts.P-1))) + 1; // the number of digits to label different bead
string format = "%0" + to_string(digits) + "d"; // format of the bead label

int ne_reject = 0, ne_num = 0;

int main() {
    // judge the calculation type
    char cal_type[16];
    constant::strtolower(Consts.calculation_type.c_str(), cal_type);

    time_t t_start = time(0); // record the start time of the program
    if (Consts.set_rc || strcmp(cal_type, "md") == 0 || strcmp(cal_type, "hmc") == 0)
        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation

    timespec tp;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp); // interval is 1e-9s
    long long seed = tp.tv_nsec; // random seed for the simulation

    // open output file
    ofstream output(Consts.out_file);
    // set digits after point in the output number
    output.setf(ios::fixed);
    output.setf(ios::showpoint);

    box B{Cell.elements, Cell.lattice_vector, Cell.N_atoms, &Cell.electron_number, seed}; // transfrom system from class cell to class box

    // set beads for path integral simulation
    if (Consts.P > 1){
        set_beads(&B); // set bead atom

        // if there's beads file, then coordinate of each bead can be set from the beginning
        if (strlen(Consts.initial_beads.c_str())){
            string initial_beads = Consts.work_path + "/" + Consts.initial_beads;
            ifstream ifbeads(initial_beads.c_str(), ios::in);
            init_beads(&B, ifbeads);
        }
    }

    // judge the potential energy and force calculator
    double (*calculate)(box *);
    char pot_type[16];
    constant::strtolower(Consts.potential_type.c_str(), pot_type);
    if (strcmp(pot_type, "dp") == 0){
        if (!strlen(Consts.DP_model_prefix.c_str())){
            canonical_dp.init(Consts.DP_model);
            dp = &canonical_dp;
	        // if the parameters of the DP potential don't contain the electron number, set the DP potential type back to the normal one
            if (!dp->dim_fparam())
                Consts.GC_type = "Discrete";
        }
        else if (strcmp(Consts.GC_type.c_str(), "Continuous") == 0){
            grand_canonical_dp.init(Consts.DP_model_prefix+".pb");
            dp = &grand_canonical_dp;
        }
        calculate = deep_potential::calculate;
    }
    else if (strcmp(pot_type, "vasp") == 0){
        // create folders for first principle calculation
        string command = ":";
        for (int i = 0 ; i < Consts.P ; ++i){
            char index[36]; // 'char *index;' will result in an error as Segmentation fault, on account of without the cooresponding memory space
            sprintf(index, format.c_str(), i);
            string folder = Consts.work_path + "_" + index;
            string command = "mkdir " + folder;
            command += " && ln -s ../" + Consts.work_path + "/INCAR " + folder + "/INCAR";
            command += " && ln -s ../" + Consts.work_path + "/KPOINTS " + folder + "/KPOINTS";
            command += " && ln -s ../" + Consts.work_path + "/POTCAR " + folder + "/POTCAR";
            system(command.c_str());
        }
        calculate = vasp::calculate;
    }
    else if (strcmp(pot_type, "abacus") == 0){
        // create folders for first principle calculation
        string command = ":";
        for (int i = 0 ; i < Consts.P ; ++i){
            char index[36]; // 'char *index;' will result in an error as Segmentation fault, on account of without the cooresponding memory space
            sprintf(index, format.c_str(), i);
            string folder = Consts.work_path + "_" + index;
            string command = "mkdir " + folder;
            command += " && ln -s ../" + Consts.work_path + "/INPUT " + folder + "/INPUT";
            command += " && ln -s ../" + Consts.work_path + "/KPT " + folder + "/KPT";
            command += " && mkdir " + folder + "/OUT.ABACUS";
            command += " && cp " + Consts.work_path + "/OUT.ABACUS/SPIN1_CHG " + folder + "/OUT.ABACUS";
            system(command.c_str());
        }
        calculate = abacus::calculate;
    }

    // model deviation settings
    model_devi *Model_Devi;
    bool perform_model_devi = Consts.model_devi_DP_models.size(); // whether to perform model deviation
    if (perform_model_devi){
        Model_Devi = new model_devi(&B);
        Model_Devi->compute_model_devi(0);
    }

    if (strcmp(cal_type, "md") == 0){
        output.precision(6);
        MD_NVE md{&B, calculate}; // MD simulation in the microcanonical ensemble
        for (int i = 0 ; i < Consts.steps ; ++i){
            // output current quantity before evolution
            output << B.kinetic_energy << " " << B.potential_energy << " " << B.energy << endl;
            md.evolve(Consts.Delta_t);

            // save current structure each certain steps
            if ((i+1) % Consts.dump_steps == 0)
                save_stru(i+1, &Cell, &B);

            // if required, perform the model deviation
            if (perform_model_devi){
                if ((i+1) % Consts.model_devi_interval == 0)
                    Model_Devi->compute_model_devi(i+1);
            }
        }
    }
    else{
        int reject = 0; // the number of rejected steps
        double ene_ave = 0; // potential energy average
        double ene_2_ave = 0; // the average of the square of potential energy
        int current_PI_step = 0; // current step of move, 0 refers to choose whether centroid move or internal move
        unordered_map<string, double> output_quantities = {{"ke", 0}, {"pe", 0}, {"etotal", 0}, {"eint", 0}, {"rc", 0}, {"mfl", 0}, 
                                                             {"mfr", 0}, {"mf", 0}, {"pe_ave", 0}, {"pe_var", 0}, {"ne", 0}};
        if (strlen(Consts.DP_model_prefix.c_str())){
            vector<string> unsupported_types{"md", "cmc"};
            for (vector<string>::const_iterator type=unsupported_types.begin() ; type < unsupported_types.end() ; ++type){
                if (strcmp(cal_type, (*type).c_str()) == 0)
                    throw runtime_error(Consts.calculation_type + " is not supported in the grand canonical simulations.");
            }
            if (strcmp(cal_type, "hmc") == 0)
            {
                if (Consts.P > 1)
                    throw invalid_argument("The path integral algorithm isn't supported in the HMC simulation.");
                if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0)
                    throw invalid_argument("Only the DP potential involving electron number is supported.");
                HMC_NVT *pmc;
                HMC_muVT mc{&B, calculate};
                pmc = &mc;
                for (int i = 0 ; i < Consts.steps ; ++i){
                    // update phyiscal quantities to new values
                    update_physical_quantities(pmc, &B, output_quantities, i, ene_ave);

                    // output current quantity before evolution
                    output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                    // evolve with certain type according to corresponding rate
                    pmc->evolve(reject);

                    // save current structure each certain steps
                    if ((i+1) % Consts.dump_steps == 0)
                        save_stru(i+1, &Cell, &B);

                    // if required, perform the model deviation
                    if (perform_model_devi){
                        if ((i+1) % Consts.model_devi_interval == 0)
                            Model_Devi->compute_model_devi(i+1);
                    }
                }
            }
            else {
                PI_HMC_muVT *pmc;
                PI_HMC_muVT mc{&B, calculate};
                if (!Consts.set_rc)
                    Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                pmc = &mc;
                if (Consts.P > 1){
                    int repeat_num = ceil((double)Consts.P/(double)Consts.J); // set continuous repeat times of the staging step
                    double target_ratio = (double)Consts.P/(double)Consts.J; // the target ratio of the staging step with relative to the centroid step
                    double staging_proportion = 1/(1+repeat_num/target_ratio); // set proportion of the staging step

                    // set proportion of the electron number change step
                    double Ne_proportion = 1/(1+((1-Consts.Ne_proportion)*(repeat_num+target_ratio))/(Consts.Ne_proportion*(target_ratio+1)*repeat_num));

                    // constrained PIHMC simulation in the canonical ensemble
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        if (strcmp(Consts.RC_type.c_str(), "Difference") == 0){
                            // update phyiscal quantities to new values
                            update_physical_quantities(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);
                        }
                        else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0){
                            // update phyiscal quantities to new values
                            update_physical_quantities_distance(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);
                        }

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (current_PI_step){
                            pmc->evolve(reject);
                            ++current_PI_step;
                            if (current_PI_step == repeat_num) current_PI_step = 0;
                        }
                        else if (B.Random() < Ne_proportion){
                            pmc->evolve(reject, true);
                            ne_num++;
                        }
                        else if (B.Random() < staging_proportion){
                            current_PI_step = 1;
                            pmc->evolve(reject);
                        }
                        else if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, false, true);
                        else
                            pmc->evolve(reject, false, true, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
                else {
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        if (strcmp(Consts.RC_type.c_str(), "Difference") == 0){
                            // update phyiscal quantities to new values
                            update_physical_quantities(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);
                        }
                        else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0){
                            // update phyiscal quantities to new values
                            update_physical_quantities_distance(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);
                        }

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (B.Random() < Consts.Ne_proportion){
                            pmc->evolve(reject, true);
                            ne_num++;
                        }
                        else if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, false, true);
                        else
                            pmc->evolve(reject, false, true, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
            }
            cout << setiosflags(ios::fixed); // set cout digits after point
            cout << "Electron Number Acceptance Rate : " << setprecision(4) << (ne_num-ne_reject)/(double)ne_num << endl;
        }
        else if (strcmp(cal_type, "hmc") == 0){
            HMC_NVT hmc(&B, calculate); // HMC simulation in the canonical ensemble
            for (int i = 0 ; i < Consts.steps ; ++i){
                // update phyiscal quantities to new values
                update_physical_quantities(&hmc, &B, output_quantities, i, ene_ave);

                // output current quantity before evolution
                output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                // evolve with the HMC scheme
                hmc.evolve(reject);

                // save current structure each certain steps
                if ((i+1) % Consts.dump_steps == 0)
                    save_stru(i+1, &Cell, &B);

                // if required, perform the model deviation
                if (perform_model_devi){
                    if ((i+1) % Consts.model_devi_interval == 0)
                        Model_Devi->compute_model_devi(i+1);
                }
            }
        }
        else if (strcmp(cal_type, "chmc") == 0){
            if (Consts.P > 1){
		        int repeat_num = ceil((double)Consts.P/(double)Consts.J); // set continuous repeat times of the staging step
		        double staging_proportion = 1/(1+repeat_num/(double)Consts.P*(double)Consts.J); // set proportion of the staging step

                // constrained PIHMC simulation in the canonical ensemble
                PI_HMC_NVT *pmc;
                if (strcmp(Consts.RC_type.c_str(), "Difference") == 0){
                    PI_HMC_NVT mc{&B, calculate};
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    pmc = &mc;
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (current_PI_step){
                            pmc->evolve(reject);
                            ++current_PI_step;
                            if (current_PI_step == repeat_num) current_PI_step = 0;
                        }
                        else if (B.Random() < staging_proportion){
                            current_PI_step = 1;
                            pmc->evolve(reject);
                        }
                        else if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, true);
                        else
                            pmc->evolve(reject, true, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
                else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0){
                    PI_HMC_NVT_DIS mc(&B, calculate);
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    pmc = &mc;
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities_distance(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (current_PI_step){
                            pmc->evolve(reject);
                            ++current_PI_step;
                            if (current_PI_step == repeat_num) current_PI_step = 0;
                        }
                        else if (B.Random() < staging_proportion){
                            current_PI_step = 1;
                            pmc->evolve(reject);
                        }
                        else if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, true);
                        else
                            pmc->evolve(reject, true, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
            }
            else {
                // constrained HMC simulation in the canonical ensemble
                HMC_NVT_RC *pmc;
                if (strcmp(Consts.RC_type.c_str(), "Difference") == 0){
                    HMC_NVT_RC mc{&B, calculate};
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    pmc = &mc;

                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, false);
                        else
                            pmc->evolve(reject, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
                else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0){
                    HMC_NVT_DIS mc(&B, calculate);
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    pmc = &mc;

                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities_distance(pmc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                        // evolve with certain type according to corresponding rate
                        if (B.Random() < Consts.hmc_proportion)
                            pmc->evolve(reject, false);
                        else
                            pmc->evolve(reject, true);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
            }
        }
        else if (strcmp(cal_type, "cmc") == 0){
            if (Consts.P > 1){
		        int repeat_num = ceil((double)Consts.P/(double)Consts.J); // set continuous repeat times of the staging step
		        double staging_proportion = 1/(1+repeat_num/(double)Consts.P*(double)Consts.J); // set proportion of the staging step
                
                PI_MC_NVT mc{&B, calculate}; // constrained PIMC simulation in the canonical ensemble
                if (!Consts.set_rc)
                    Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                for (int i = 0 ; i < Consts.steps ; ++i){
                    // update phyiscal quantities to new values
                    update_physical_quantities(&mc, &B, output_quantities, i, ene_ave, ene_2_ave);

                    // output current quantity before evolution
                    output_physical_quantities(output, Consts.quantity_digits, output_quantities);

                    // evolve with certain type according to corresponding rate
                    if (current_PI_step){
                        mc.evolve(reject);
                        ++current_PI_step;
                        if (current_PI_step == repeat_num) current_PI_step = 0;
                    }
                    else if (B.Random() < staging_proportion){
                        current_PI_step = 1;
                        mc.evolve(reject);
                    }
                    else
                        mc.evolve(reject, true);

                    // save current structure each certain steps
                    if ((i+1) % Consts.dump_steps == 0)
                        save_stru(i+1, &Cell, &B);

                    // if required, perform the model deviation
                    if (perform_model_devi){
                        if ((i+1) % Consts.model_devi_interval == 0)
                            Model_Devi->compute_model_devi(i+1);
                    }
                }
            }
            else {
                if (strcmp(Consts.RC_type.c_str(), "Difference") == 0){
                    MC_NVT_RC mc{&B, calculate}; // constrained MC simulation in the canonical ensemble
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities(&mc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);
                        mc.evolve(reject);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
                else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0){
                    MC_NVT_DIS mc{&B, calculate}; // constrained MC simulation in the canonical ensemble
                    if (!Consts.set_rc)
                        Consts.Print("ALL_INPUT"); // output values of all parameters used in the simulation
                    for (int i = 0 ; i < Consts.steps ; ++i){
                        // update phyiscal quantities to new values
                        update_physical_quantities_distance(&mc, &B, output_quantities, i, ene_ave, ene_2_ave);

                        // output current quantity before evolution
                        output_physical_quantities(output, Consts.quantity_digits, output_quantities);
                        mc.evolve(reject);

                        // save current structure each certain steps
                        if ((i+1) % Consts.dump_steps == 0)
                            save_stru(i+1, &Cell, &B);

                        // if required, perform the model deviation
                        if (perform_model_devi){
                            if ((i+1) % Consts.model_devi_interval == 0)
                                Model_Devi->compute_model_devi(i+1);
                        }
                    }
                }
            }
        }
        cout << setiosflags(ios::fixed); // set cout digits after point
        cout << "Acceptance Rate : " << setprecision(4) << (Consts.steps-reject)/(double)Consts.steps << endl;
        cout << "Average Potential Energy : " << setprecision(6) << ene_ave << endl;
    }
    B.regular();
    Cell.print_stru_file("STRU_after"); // save last structure

    // remove folders for first principle calcuation
    if (strcmp(pot_type, "vasp") == 0 || strcmp(pot_type, "abacus") == 0){
        string command = "rm " + Consts.work_path + "_* -r";
        system(command.c_str());
    }

    time_t t_end = time(0); // record the end time of the program
    run_time(t_start, t_end); // print running time of the program
    if (perform_model_devi)
        delete Model_Devi;
    return 0;
}
