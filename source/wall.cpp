#include "wall.h"
using namespace cpihmc;
using namespace std;

cpihmc::wall::wall(box *Box, const unordered_map<index_t, prec_t> &Walls):hard_bdry(Box), WallNum(Walls.size()), WallJudges({&wall::wall_xlo, &wall::wall_xhi, &wall::wall_ylo, &wall::wall_yhi, &wall::wall_zlo, &wall::wall_zhi}), WallPoses(6)
{
    if (WallNum)
        for (const auto &Wall : Walls)
        {
            WallPoses[Wall.first] = Wall.second; // set the position(s) of the set wall(s)
            UsedJudges.push_back(WallJudges[Wall.first]); // set the judge criterion(s) of the set wall(s)
        }
}

const bool_t cpihmc::wall::judge_accept() const
{
    if (WallNum)
    {
        for (const auto &Index : Box->AtomIndexes) for (const auto &Judge : UsedJudges) if ((this->*Judge)(Box->Atoms[Index].Coord)) return false;
        if (Box->NBead > 1) for (const auto &Bead : Box->Beads) for (index_t i = 0 ; i < Box->NBead ; ++i) for (const auto &Judge : UsedJudges) if ((this->*Judge)(Bead.Coords[i])) return false;
    }
    return true;
}