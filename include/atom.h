#ifndef CPIHMC_ATOM_H
#define CPIHMC_ATOM_H

#include <vector>
#include <stdexcept>
#include "type.h"
#include "vec3.h"

namespace cpihmc
{
    class atom
    {
        public:
            vec3<prec_t> Coord;
            vec3<prec_t> Vel;
            vec3<prec_t> Force;
            vec3<bool_t> Move;
            const prec_t Mass;
        public:
            explicit atom(const prec_t Mass):Vel(0.0, 0.0, 0.0), Force(0.0, 0.0, 0.0), Move(true, true, true), Mass(Mass){}
            explicit atom(const vec3<prec_t> &Coord, const prec_t Mass):Coord(Coord), Vel(0.0, 0.0, 0.0), Force(0.0, 0.0, 0.0), Move(true, true, true), Mass(Mass){}
            explicit atom(const vec3<prec_t> &Coord, const vec3<bool_t> &Move, const prec_t Mass):Coord(Coord), Vel(0.0, 0.0, 0.0), Force(0.0, 0.0, 0.0), Move(Move), Mass(Mass){}
            inline const atom &operator=(const atom &Atom);
            ~atom() = default;
    };

    class bead
    {
        public:
            atom * const BeadAtom;
            std::vector<vec3<prec_t>> Coords;
            prec_t KinEng;
            const prec_t Mass;
        public:
            bead(atom * const Atom, const index_t P, const prec_t Mass):BeadAtom(Atom), Coords(P, BeadAtom->Coord), KinEng(0.0), Mass(Mass){}
            ~bead() = default;
            void back_centroid();
    };
}

const cpihmc::atom &cpihmc::atom::operator=(const cpihmc::atom &Atom)
{
    if (this != &Atom)
    {
        if (Mass != Atom.Mass) throw std::invalid_argument("The properties of an atom cannot be assigned to another atom with a different mass.");
        Coord = Atom.Coord;
        Vel = Atom.Vel;
        Force = Atom.Force;
        Move = Atom.Move;
    }
    return *this;
}

#endif