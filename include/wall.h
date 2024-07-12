//
// Created by Jin Bin on 2023/04/17.
// to set wall in the the cell
//

#ifndef CPIHMC_WALL_H
#define CPIHMC_WALL_H

#include "box.h"

class hard_boundary // base class
{
public:
    hard_boundary(box *box):Box(box){}
    virtual bool judge_accept() const = 0;
    virtual ~hard_boundary() = default;
protected:
    box *Box; // the simulated system
};

class wall;
typedef bool (wall::*wall_type)(const atom & atom) const; // define a member function pointer type

class wall : public hard_boundary {
public:
    wall(box *);
    // judge whether to accept the structure based on the current atom coordinates relative to the set wall(s)
    bool judge_accept() const;
    ~wall() = default;
private:
    int n_wall = Consts.wall.size(); // the number of the set wall(s)
    vector<double> wall_positions = vector<double>(6); // the positions of the wall(s)
    bool wall_check = true; // whether the wall check is needed
    static map<string, int> wall_index; // the indexes of different wall types
    // the judge criterions of different wall types
    map<string, wall_type> wall_judge = {{"xlo", &wall::wall_xlo}, {"xhi", &wall::wall_xhi}, 
    {"ylo", &wall::wall_ylo}, {"yhi", &wall::wall_yhi}, {"zlo", &wall::wall_zlo}, {"zhi", &wall::wall_zhi}};
    vector<wall_type> used_judge; // the judge criterions of the set wall(s)
    bool wall_xlo(const atom &atom) const
    // whether the atom exceeds the x low wall
    {
        return atom.r.x < wall_positions[0];
    }
    bool wall_xhi(const atom &atom) const
    // whether the atom exceeds the x high wall
    {
        return atom.r.x > wall_positions[1];
    }
    bool wall_ylo(const atom &atom) const
    // whether the atom exceeds the y low wall
    {
        return atom.r.y < wall_positions[2];
    }
    bool wall_yhi(const atom &atom) const
    // whether the atom exceeds the y high wall
    {
        return atom.r.y > wall_positions[3];
    }
    bool wall_zlo(const atom &atom) const
    // whether the atom exceeds the z low wall
    {
        return atom.r.z < wall_positions[4];
    }
    bool wall_zhi(const atom &atom) const
    // whether the atom exceeds the z high wall
    {
        return atom.r.z > wall_positions[5];
    }
};

#endif //CPIHMC_WALL_H
