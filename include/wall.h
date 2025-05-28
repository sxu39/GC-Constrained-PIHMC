#ifndef CPIHMC_WALL_H
#define CPIHMC_WALL_H

#include "hard_bdry.h"

namespace cpihmc
{
    class wall;
    typedef const bool_t (wall::*wall_type)(const vec3<prec_t> &) const; // define a member function pointer type

    class wall : public hard_bdry
    {
        private:
            const size_t WallNum; // the number of the set wall(s)
            const std::vector<wall_type> WallJudges; // the judge criterions of different wall types
            std::vector<wall_type> UsedJudges; // the judge criterion(s) of the set wall(s)
            std::vector<prec_t> WallPoses; // the positions of the wall(s)
        private:
            const bool_t wall_xlo(const vec3<prec_t> &Coord) const {return Coord.x < WallPoses[0];} // whether the atom exceeds the x low wall
            const bool_t wall_xhi(const vec3<prec_t> &Coord) const {return Coord.x > WallPoses[1];} // whether the atom exceeds the x high wall
            const bool_t wall_ylo(const vec3<prec_t> &Coord) const {return Coord.y < WallPoses[2];} // whether the atom exceeds the y low wall
            const bool_t wall_yhi(const vec3<prec_t> &Coord) const {return Coord.y > WallPoses[3];} // whether the atom exceeds the y high wall
            const bool_t wall_zlo(const vec3<prec_t> &Coord) const {return Coord.z < WallPoses[4];} // whether the atom exceeds the z low wall
            const bool_t wall_zhi(const vec3<prec_t> &Coord) const {return Coord.z > WallPoses[5];} // whether the atom exceeds the z high wall
        public:
            explicit wall(box *, const std::unordered_map<index_t, prec_t> &);
            ~wall() = default;
            const bool_t judge_accept() const final; // judge whether to accept the structure based on the current atom coordinates relative to the set wall(s)
    };
}

#endif