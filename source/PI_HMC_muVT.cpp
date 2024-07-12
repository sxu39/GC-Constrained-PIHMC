#include "PI_HMC_muVT.h"
#include "gen_func.h"

PI_HMC_muVT::PI_HMC_muVT(box *Box, double (*calculate)(box *)){
    this->Box = Box;
    this->calculate = calculate;
    system("mkdir ALL_STRU");
    if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0){
        init_models();
        set_model();
    }
    // calculate system energy and force for the first time
    Box->potential_energy = calculate(Box);
    Box->energy = Box->kinetic_energy + Box->potential_energy;
    if (strcmp(Consts.RC_type.c_str(), "Difference") == 0)
        PI_HMC_NVT::set_constraint_atoms();
    else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0)
        set_constraint_atoms();
    if (!Consts.set_rc){
        if (strcmp(Consts.RC_type.c_str(), "Difference") == 0)
            Consts.reaction_coordinate = ::react_coord(hidari_atom, naka_atom, migi_atom);
        else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0)
            Consts.reaction_coordinate = ::react_coord(correlation_atom, constraint_atom);
    }
    react_coord = Consts.reaction_coordinate;
    hardBoundary = new wall {Box};
}

void PI_HMC_muVT::evolve(int &current, bool Ne, bool centroid, bool MC_step) const{
    if (Ne)
        electron_number_evolve(current);
    else if (centroid)
        centroid_evolve(current, MC_step);
    else
        internal_evolve(current);
}

void PI_HMC_muVT::init_models()
{
    string model_name;
    double Ne;
    int model_number = round((Ne_range.second - Ne_range.first)/delta_Ne) + 1;
    dp_models.resize(model_number);
    for (int i = 0 ; i < model_number ; ++i){
        Ne = Ne_range.first+delta_Ne*i;
        stringstream sNe;
        sNe << setiosflags(ios::fixed) << setprecision(2) << Ne;
        model_name = Consts.DP_model_prefix + sNe.str() + ".pb";
        dp_models[i].init(model_name);
    }
}

void PI_HMC_muVT::HMC_evolve(int &current) const {
    if (strcmp(Consts.RC_type.c_str(), "Difference") == 0)
        PI_HMC_NVT::HMC_evolve(current);
    else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0)
        PI_HMC_NVT_DIS::HMC_evolve(current);
}

void PI_HMC_muVT::MC_evolve(int &current) const {
    if (strcmp(Consts.RC_type.c_str(), "Difference") == 0)
        PI_HMC_NVT::MC_evolve(current);
    else if (strcmp(Consts.RC_type.c_str(), "Distance") == 0)
        PI_HMC_NVT_DIS::MC_evolve(current);
}

void PI_HMC_muVT::electron_number_evolve(int &current) const {
    double temp_energy = Box->potential_energy;
    double temp_electron_number = *Box->electron_number;
    double factor = 1;
    if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0){
        if (Box->Random() < 0.5)
            *Box->electron_number = min(*Box->electron_number+delta_Ne, Ne_range.second);
        else
            *Box->electron_number = max(*Box->electron_number-delta_Ne, Ne_range.first);
    }
    else if (strcmp(Consts.GC_type.c_str(), "Continuous") == 0){
        pair<double, double> current_range, trial_range;
        current_range.first = max(*Box->electron_number-0.5*delta_Ne, Ne_range.first);
        current_range.second = min(*Box->electron_number+0.5*delta_Ne, Ne_range.second);
        *Box->electron_number = current_range.first + Box->Random()*(current_range.second-current_range.first);
        trial_range.first = max(*Box->electron_number-0.5*delta_Ne, Ne_range.first);
        trial_range.second = min(*Box->electron_number+0.5*delta_Ne, Ne_range.second);
        factor = (current_range.second-current_range.first) / (trial_range.second-trial_range.first);
    }
    if (*Box->electron_number != temp_electron_number){
        if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0)
            set_model();
        Box->potential_energy = calculate(Box);
        if (Box->Random() > min(1.0, factor*exp(((*Box->electron_number-temp_electron_number)*mu-(Box->potential_energy-temp_energy))/(Box->k_B*Box->T)))){
            *Box->electron_number = temp_electron_number;
            if (strcmp(Consts.GC_type.c_str(), "Discrete") == 0)
                set_model();
            Box->potential_energy = calculate(Box);
            current++;
            ne_reject++;
        }
        Box->energy = Box->kinetic_energy + Box->potential_energy;
    }
}

void PI_HMC_muVT::set_model() const {
    int index = round((*Box->electron_number-Ne_range.first)/delta_Ne);
    dp = const_cast<deepmd::DeepPot *>(&dp_models[index]);
}
