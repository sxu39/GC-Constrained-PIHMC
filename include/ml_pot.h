#ifndef CPIHMC_ML_POT_H
#define CPIHMC_ML_POT_H

#include "pot.h"
#include "box.h"
#include "deepmd/DeepPot.h"

namespace cpihmc
{
    class deep_pot : public pot
    {
        friend class model_devi;
        private:
            deepmd::DeepPot DeepPot;
        public:
            explicit deep_pot(const std::string &DeepPotName){DeepPot.init(DeepPotName);}
            ~deep_pot() = default;
            prec_t infer(box *) final;
        private:
            static void provide_structures(std::vector<prec_t> &, std::vector<index_t> &, std::vector<prec_t> &, box *);
            void single_point_calculate(std::vector<prec_t> &, std::vector<prec_t> &, const std::vector<prec_t> &, const std::vector<index_t> &, const std::vector<prec_t> &);
            void single_point_calculate(std::vector<prec_t> &, std::vector<prec_t> &, const std::vector<prec_t> &, const std::vector<index_t> &, const std::vector<prec_t> &, const std::vector<prec_t> &);
    };
}

#endif
