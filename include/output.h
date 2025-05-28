#ifndef CPIHMC_OUTPUT_H
#define CPIHMC_OUTPUT_H

#include "box.h"
#include "rxn_coord.h"
#include <unordered_map>

namespace cpihmc
{
    template <class T> typename std::enable_if<std::is_same<T, std::string>::value, const char *>::type filter(const T &Str){return Str.c_str();}
    template <class T> typename std::enable_if<!std::is_same<T, std::string>::value, const T &>::type filter(const T &Str){return Str;}
    template <class ... T> const std::string format(const char *Fmt, const T & ... Args);

    class output
    {
        private:
            const input &Input;
            const box &Box;
            std::ofstream OutFile;
            std::unordered_map<std::string, prec_t> OutQuant;
            size_t StepColWidth;
            size_t ColWidth;
            static std::unordered_map<std::string, std::string> PhyQuantName;
        private:
            void insert_rxn_coord_quant(const std::string &Suffix="");
            const prec_t internal_energy_estimator() const;
            void update_rxn_coord_quant(const rxn_coord * const, const std::string &Suffix="");
            void print_stru_file(const std::string &, const index_t &NSpin=1, const bool_t &Direct=false, const bool_t &Vel=false,
                                 const bool_t &Magmom=false) const;
            void print_beads_file(const std::string &) const;
        public:
            size_t RejectTimes;
            time_t BeginTime, EndTime;
        public:
            explicit output(const input &, const box &);
            ~output() = default;
            void update_physical_quantities(const index_t, const std::vector<rxn_coord *> &);
            void print_physical_quantities(const index_t);
            void save_structure(const index_t) const;
            void simulation_summary() const;
    };
}

template <class ... T> 
const std::string cpihmc::format(const char *Fmt, const T & ... Args)
{
    const index_t Size = snprintf(nullptr, 0, Fmt, filter(Args)...) + 1;
    std::string Dst(Size, '\0');
    snprintf(&Dst[0], Size, Fmt, filter(Args)...);
    Dst.pop_back();
    return Dst;
}

#endif