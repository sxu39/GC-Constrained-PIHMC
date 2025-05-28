#ifndef CPIHMC_MODEL_DEVI_H
#define CPIHMC_MODEL_DEVI_H

#include "box.h"
#include "ml_pot.h"

namespace cpihmc
{
    class model_devi
    {
        private:
            box * const Box;
            size_t ModelNum;
            std::vector<deepmd::DeepPot> Models;
            std::ofstream OutFile;
        private:
            void compute_force_virials(std::vector<std::vector<prec_t>> &, std::vector<std::vector<prec_t>> &, const std::vector<prec_t> &, const std::vector<index_t> &, const std::vector<prec_t> &);
            template <class T> void compute_avg(std::vector<T> &, const std::vector<std::vector<T>> &) const;
            template <class T> void compute_std(std::vector<T> &, const std::vector<T> &, const std::vector<std::vector<T>> &, const index_t) const;
            static void ana_st(prec_t &, prec_t &, prec_t &, const std::vector<prec_t> &);
            void cope_with_force(prec_t &, prec_t &, prec_t &, const std::vector<std::vector<prec_t>> &) const;
            void cope_with_virial(prec_t &, prec_t &, prec_t &, const std::vector<std::vector<prec_t>> &) const;
        public:
            explicit model_devi(box * const, const std::vector<std::string> &, const std::string &);
            void compute_model_devi(const index_t);
            ~model_devi();
    };
}

#endif