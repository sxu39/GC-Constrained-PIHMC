#include "HMC_muVT.h"

HMC_muVT::HMC_muVT(box *Box, double (*calculate)(box *)):HMC_NVT(Box, calculate){}

void HMC_muVT::evolve(int &current) const
{
    // evolve with certain type according to corresponding rate
    if (Box->Random() < Consts.Ne_proportion){
        electron_number_evolve(current);
        ne_num++;
    }
    else
        HMC_NVT::evolve(current);
}

void HMC_muVT::electron_number_evolve(int &current) const {
    double temp_energy = Box->potential_energy;
    double temp_electron_number = *Box->electron_number;
    double factor = 1;
    pair<double, double> current_range, trial_range;
    current_range.first = max(*Box->electron_number-0.5*delta_Ne, Ne_range.first);
    current_range.second = min(*Box->electron_number+0.5*delta_Ne, Ne_range.second);
    *Box->electron_number = current_range.first + Box->Random()*(current_range.second-current_range.first);
    trial_range.first = max(*Box->electron_number-0.5*delta_Ne, Ne_range.first);
    trial_range.second = min(*Box->electron_number+0.5*delta_Ne, Ne_range.second);
    factor = (current_range.second-current_range.first) / (trial_range.second-trial_range.first);
    if (*Box->electron_number != temp_electron_number){
        Box->potential_energy = calculate(Box);
        if (Box->Random() > min(1.0, factor*exp(((*Box->electron_number-temp_electron_number)*mu-(Box->potential_energy-temp_energy))/(Box->k_B*Box->T)))){
            *Box->electron_number = temp_electron_number;
            Box->potential_energy = calculate(Box);
            current++;
            ne_reject++;
        }
        Box->energy = Box->kinetic_energy + Box->potential_energy;
    }
}