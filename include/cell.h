#ifndef CPIHMC_CELL_H
#define CPIHMC_CELL_H

#include <cmath>
#include <sstream>
#include "const.h"
#include "input.h"
#include "element.h"
#include "mat3.h"

namespace cpihmc
{
    class cell
    {
        friend class box;
        friend class output;
        private:
            std::string *AtomLabel;
            prec_t *AtomMass;
            std::string *PseudoFileName;
            std::string *OrbitalFileName;
            prec_t *Magnetization;
            prec_t **AtomMagnet;
            vec3<prec_t> **AtomMagnetVec;
            bool_t SetElecNum;
            prec_t ElecNum;
            prec_t LatConst;
            prec_t LatConstAng;
            mat3<prec_t> LatVec;
            prec_t Omega;
            std::string CoordType;
            size_t NType;
            size_t NAtoms;
            std::vector<atom> Atoms;
            std::vector<element> Elements;
            std::vector<index_t> Masks;
        private:
            void setup_cell(const std::string &);
	        void read_atom_species(std::ifstream &); // read in the atom information for each type of atom
	        const bool_t read_atom_positions(std::ifstream &); // read in atomic positions
        public:
            explicit cell(const input &);
            ~cell();
    };
}

#endif