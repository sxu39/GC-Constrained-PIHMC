#ifndef CPIHMC_ELEMENT_H
#define CPIHMC_ELEMENT_H

#include "atom.h"

namespace cpihmc
{
    class element
    {
        private:
            std::vector<atom> &Atoms;
            size_t Begin;
        public:
            const std::string Label;
            size_t Number;
            const prec_t Mass;
        public:
            explicit element(const std::string &Label, const prec_t Mass, std::vector<atom> &Atoms):Label(Label), Number(0), Mass(Mass), Atoms(Atoms), Begin(0){}
            ~element() = default;
            inline void set_begin(const size_t Begin);
            const std::vector<atom>::const_iterator begin() const {return Atoms.cbegin()+Begin;}
            const std::vector<atom>::const_iterator end() const {return Atoms.cbegin()+Begin+Number;}
            const std::vector<atom>::iterator begin(){return Atoms.begin()+Begin;}
            const std::vector<atom>::iterator end(){return Atoms.begin()+Begin+Number;}
            inline const element &operator=(const element &);
            const atom &operator[](index_t Index) const {return Atoms[Begin+Index];}
            atom &operator[](index_t Index){return Atoms[Begin+Index];}
    };
}

void cpihmc::element::set_begin(const size_t Begin)
{
    this->Begin = Begin;
    if (this->Begin + Number > Atoms.size()) throw std::invalid_argument("There aren't enough atoms.");
}

const cpihmc::element &cpihmc::element::operator=(const cpihmc::element &Element)
{
    if (this != &Element)
    {
        if (Label != Element.Label) throw std::invalid_argument("The properties of an element cannot be assigned to another element with a different label.");
        if (Mass != Element.Mass) throw std::invalid_argument("The properties of an element cannot be assigned to another element with a different mass.");
        Atoms = Element.Atoms;
        Begin = Element.Begin;
        Number = Element.Number;
    }
    return *this;
}

#endif