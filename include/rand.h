#ifndef CPIHMC_RAND_H
#define CPIHMC_RAND_H

#include <random>
#include "type.h"

namespace cpihmc
{
    class rand
    {
        private:
            std::random_device rand_device;
            std::mt19937 rand_engine;
            std::uniform_real_distribution<prec_t> uni_dis = std::uniform_real_distribution<prec_t>(0.0, 1.0);
            std::normal_distribution<prec_t> nor_dis = std::normal_distribution<prec_t>(0.0, 1.0);
        private:
            rand():rand_engine(rand_device()){}
            ~rand() = default;
            rand(const rand &) = delete;
            rand(rand &&) = delete;
            const rand &operator=(const rand &) = delete;
        public:
            static rand &get_instance()
            {
                static rand instance;
                return instance;
            }
            prec_t uniform()
            {
                return uni_dis(rand_engine);
            }
            index_t uniform_int(index_t end)
            {
                return floor(uniform()*end);
            }
            index_t uniform_int(index_t start, index_t end)
            {
                return floor(uniform()*(end-start)+start);
            }
            prec_t normal()
            {
                return nor_dis(rand_engine);
            }
    };
}

#endif