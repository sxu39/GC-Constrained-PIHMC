#ifndef CPIHMC_INPUT_H
#define CPIHMC_INPUT_H

#include "type.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <array>
#include <unordered_map>
#include <cassert>
#include <iomanip>

namespace cpihmc
{
    bool_t scan_begin(std::ifstream &InFile, const std::string &TargetName, const bool_t Restart=true);
    template <class T> void read_value(std::ifstream &InFile, T &Val);
    template <class T> void read_value(std::ifstream &InFile, std::vector<T> &ValVec);
    inline void str2lower(const char *StrIn, char *StrOut);
    template <class T> void out_param(std::ofstream &OutFile, const std::string &Name, const T &Val, const std::string &Explan="");

    const std::unordered_map<std::string, size_t> RxnCoordAtomNums = {{"DIST", 2}, {"DIFF", 3}};
    const std::unordered_map<std::string, size_t> RxnCoordParamNums = {{"DIST", 1}, {"DIFF", 3}};
    const std::unordered_map<std::string, index_t> WallIndexes = {{"xlo", 0}, {"xhi", 1}, {"ylo", 2}, {"yhi", 3}, {"zlo", 4}, {"zhi", 5}};

    struct rxn_coord_info
    {
        std::string type;
        std::vector<index_t> atom_indexes;
        prec_t value;
        std::vector<prec_t> params;
        bool_t set_rxn_coord;
        rxn_coord_info(const std::string &type, const std::vector<index_t> &atom_indexes, const prec_t value, const std::vector<prec_t> &params, const bool_t set_rxn_coord=true):type(type), atom_indexes(atom_indexes), value(value), params(params), set_rxn_coord(set_rxn_coord){}
    };

    class input
    {
        private:
            // General
            std::string stru_file;
            std::string pot_type;
            std::string deep_pot_model;
            std::string simu_type;
            size_t steps;
            // System
            prec_t temp;
            size_t n_type;
            std::vector<std::string> element_type;
            std::vector<index_t> element_index;
            std::unordered_map<std::string, prec_t> wall;
            std::unordered_map<std::string, std::vector<index_t>> group_indexes;
            // Hybrid Monte Carlo
            size_t n_evol_step;
            prec_t time_step;
            prec_t mass_scal;
            // Constraint
            prec_t hybrid_monte_carlo_ratio;
            std::unordered_map<std::string, index_t> virt_atom;
            std::vector<std::string> virt_atom_names;
            std::vector<std::vector<index_t>> virt_atom_indexes;
            std::vector<rxn_coord_info> rxn_coord;
            // Grand Canonical Ensemble
            prec_t elec_num_ratio;
            prec_t mu;
            std::pair<prec_t, prec_t> elec_num_range;
            prec_t elec_num_width;
            // Path Integral
            std::string beads_file;
            size_t n_bead;
            size_t n_change_bead;
            std::vector<index_t> bead_index;
            // Model Deviation
            std::vector<std::string> model_devi_deep_pot_models;
            size_t model_devi_intvl;
            std::string model_devi_file;
            // Output
            std::string phy_quant_file;
            size_t col_width;
            std::vector<std::string> out_phy_quant;
            std::unordered_map<std::string, index_t> phy_quant_digits;
            size_t phy_quant_intvl;
            size_t stru_intvl;
        private:
            void insert_rxn_coord_digits(const std::string &Suffix="");
        public:
            explicit input(const std::string &);
            ~input() = default;
            void read(const std::string &);
            void print(const std::string &) const;
            inline const std::string &get_stru_file() const {return stru_file;}
            inline const std::string &get_pot_type() const {return pot_type;}
            inline const std::string &get_deep_pot_model() const {return deep_pot_model;}
            inline const std::string &get_simu_type() const {return simu_type;}
            inline const size_t get_steps() const {return steps;}
            inline const prec_t get_temp() const {return temp;}
            inline const size_t get_n_type() const {return n_type;}
            inline const std::unordered_map<std::string, index_t> get_element_type() const;
            inline const std::unordered_map<index_t, prec_t> get_wall() const;
            inline const size_t get_n_evol_step() const {return n_evol_step;}
            inline const prec_t get_time_step() const {return time_step;}
            inline const prec_t get_mass_scal() const {return mass_scal;}
            inline const prec_t get_hybrid_monte_carlo_ratio() const {return hybrid_monte_carlo_ratio;}
            inline const std::vector<std::vector<index_t>> &get_virt_atom_indexes() const {return virt_atom_indexes;}
            inline const std::vector<rxn_coord_info> &get_rxn_coord() const {return rxn_coord;}
            inline const prec_t get_elec_num_ratio() const {return elec_num_ratio;}
            inline const prec_t get_mu() const {return mu;}
            inline const std::pair<prec_t, prec_t> &get_elec_num_range() const {return elec_num_range;}
            inline const prec_t get_elec_num_width() const {return elec_num_width;}
            inline const std::string &get_beads_file() const {return beads_file;}
            inline const size_t get_n_bead() const {return n_bead;}
            inline const size_t get_n_change_bead() const {return n_change_bead;}
            inline const std::vector<index_t> &get_bead_index() const {return bead_index;}
            inline const std::vector<std::string> &get_model_devi_deep_pot_models() const {return model_devi_deep_pot_models;}
            inline const size_t get_model_devi_intvl() const {return model_devi_intvl;}
            inline const std::string &get_model_devi_file() const {return model_devi_file;}
            inline const std::string &get_phy_quant_file() const {return phy_quant_file;}
            inline const size_t get_col_width() const {return col_width;}
            inline const std::vector<std::string> &get_out_phy_quant() const {return out_phy_quant;}
            inline const std::unordered_map<std::string, index_t> &get_phy_quant_digits() const {return phy_quant_digits;}
            inline const size_t get_phy_quant_intvl() const {return phy_quant_intvl;}
            inline const size_t get_stru_intvl() const {return stru_intvl;}
    };
}

template <class T>
void cpihmc::read_value(std::ifstream &InFile, T &Val)
{
    InFile >> Val;
    std::string Line;
    getline(InFile, Line);
}

template <class T>
void cpihmc::read_value(std::ifstream &InFile, std::vector<T> &ValVec)
{
    T Temp;
    while (InFile.good() && InFile.peek() != '\n')
    {
        if (InFile.peek() == ' ' || InFile.peek() == '\t') InFile.get();
        else if (InFile.peek() == '#'){InFile.ignore(150, '\n'); break;}
        else {InFile >> Temp; ValVec.push_back(Temp);}
    }
    if (InFile.peek() == '\n') InFile.get(); // get rid of '\n'
}

void cpihmc::str2lower(const char *StrIn, char *StrOut)
{
    size_t Len = strlen(StrIn);
    for (index_t i = 0 ; i < Len ; ++i)
        StrOut[i] = tolower(StrIn[i]);
    StrOut[Len] = '\0';
}

template <class T>
void cpihmc::out_param(std::ofstream &OutFile, const std::string &Name, const T &Val, const std::string &Explan){OutFile << std::setw(32) << Name << Val << " # " << Explan << std::endl;}

const std::unordered_map<std::string, cpihmc::index_t> cpihmc::input::get_element_type() const
{
    assert(n_type == element_type.size());
    std::unordered_map<std::string, index_t> ElementType;
    if (element_index.size())
        for (index_t i = 0 ; i < n_type ; ++i)
            ElementType.insert(std::pair<std::string, index_t>(element_type[i], element_index[i]));
    else
        for (index_t i = 0 ; i < n_type ; ++i)
            ElementType.insert(std::pair<std::string, index_t>(element_type[i], i));
    return ElementType;
}

const std::unordered_map<cpihmc::index_t, cpihmc::prec_t> cpihmc::input::get_wall() const
{
    std::unordered_map<index_t, prec_t> Walls;
    for (const auto &Wall : wall) Walls.insert({WallIndexes.at(Wall.first), Wall.second});
    return Walls;
}

#endif