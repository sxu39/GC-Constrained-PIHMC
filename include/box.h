#ifndef CPIHMC_BOX_H
#define CPIHMC_BOX_H

#include "cell.h"
#include "rand.h"

namespace cpihmc
{
    class box
    {
        friend class output;
        private:
            cell Cell;
        private:
            static const index_t find_ge_index(const std::vector<index_t> &, const index_t);
            const index_t element2atom(const index_t, const index_t) const;
            void update_kinetic_energy();
        public:
            mat3<prec_t> LatVec;
            prec_t &ElecNum;
            const prec_t MassScal;
            const size_t NType;
            const std::unordered_map<std::string, index_t> ElementType;
            const size_t &NAtoms;
            std::vector<atom> &Atoms;
            std::vector<element> &Elements;
            std::vector<index_t> &Masks;
            std::vector<index_t> AtomIndexes;
            const size_t NBead;
            const std::vector<index_t> BeadIndexes;
            std::vector<bead> Beads;
            const prec_t Temp;
            prec_t KinEng;
            prec_t PotEng;
            prec_t TotEng;
            prec_t QtmKinEng;
        public:
            explicit box(const input &Input);
            ~box() = default;
            void update_total_energy(const bool_t UpdateKinEng=false){if (UpdateKinEng) update_kinetic_energy(); TotEng = KinEng + PotEng;}
            void insert_atom(const index_t);
            void insert_atom(const index_t, const index_t);
            void remove_atom(const index_t);
            void remove_atom(const index_t, const index_t);
    };
}

#endif