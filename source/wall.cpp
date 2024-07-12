//
// Created by Jin Bin on 2023/04/17.
//

#include "wall.h"

map<string, int> wall::wall_index = {{"xlo", 0}, {"xhi",1}, {"ylo",2}, {"yhi",3}, {"zlo",4}, {"zhi",5}};

wall::wall(box *box):hard_boundary(box){
    if (n_wall){
        for (int i = 0 ; i < n_wall ; ++i){
            // set the position(s) of the set wall(s)
            wall_positions[wall_index[Consts.wall[i]]] = Consts.wall_pos[i];

            // set the judge criterion(s) of the set wall(s)
            used_judge.push_back(wall_judge[Consts.wall[i]]);
        }
    }
    else
        wall_check = false;
}

bool wall::judge_accept() const{
    if (n_wall)
        for (int i = 0 ; i < Box->ntype ; ++i){
            for (int j = 0 ; j < Box->elements[i].number ; ++j){
                for (vector<wall_type>::const_iterator wj = used_judge.begin() ; wj < used_judge.end() ; ++wj){
                    if ((this->**wj)(Box->elements[i].atoms[j]))
                        return false;
                }
            }
        }
    return true;
}