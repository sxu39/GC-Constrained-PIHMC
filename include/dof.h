#ifndef CPIHMC_DOF_H
#define CPIHMC_DOF_H

#include "atom.h"
#include <numeric>

namespace cpihmc
{
    class dof_vec3 : public vec3<prec_t>
    {
        private:
            size_t VecNum;
            std::vector<vec3<prec_t> *> AllVecs;
        public:
            explicit dof_vec3(const std::vector<vec3<prec_t> *> &AllVecs):VecNum(AllVecs.size()), AllVecs(AllVecs){synchronize();}
            const dof_vec3 &operator=(const vec3<prec_t> &u) final
            {
                vec3<prec_t> DiffVec = u - vec3<prec_t>{x, y, z};
                x = u.x; y = u.y; z = u.z;
                for (auto &Vec : AllVecs) *Vec += DiffVec;
                return *this;
            }
            dof_vec3 &operator+=(const vec3<prec_t> &u) final {vec3<prec_t>::operator+=(u); for (auto &Vec : AllVecs) *Vec += u; return *this;}
            dof_vec3 &operator-=(const vec3<prec_t> &u) final {vec3<prec_t>::operator-=(u); for (auto &Vec : AllVecs) *Vec -= u; return *this;}
            void synchronize()
            {
                vec3<prec_t> VecSum;
                for (auto &Vec : AllVecs) VecSum += *Vec;
                vec3<prec_t>::operator=(VecSum / (prec_t)VecNum);
            }
    };

    const std::vector<vec3<prec_t> *> obtain_coord_pointer(const std::vector<atom *> &AllAtoms);
    const prec_t obtain_mean_mass(const std::vector<atom *> &AllAtoms);

    class dof_atom
    {
        friend class dof;
        friend class virt_atom;
        private:
            std::vector<atom *> AllAtoms;
            dof_vec3 Coord;
            vec3<prec_t> Vel;
            vec3<prec_t> Force;
            const prec_t Mass;
        private:
            void update_force();
        public:
            explicit dof_atom(const std::vector<atom *> &AllAtoms);
            ~dof_atom() = default;
    };

    class dof
    {
        private:
            dof_atom * const VirtAtom;
            vec3<prec_t> &Force;
        public:
            vec3<prec_t> &Coord;
            vec3<prec_t> &Vel;
            const prec_t &Mass;
            const size_t AtomNum;
        public:
            explicit dof(atom &Atom);
            explicit dof(dof_atom &VirtAtom);
            ~dof() = default;
            const vec3<prec_t> get_force() const;
            void synchronize() const;
    };
}

#endif