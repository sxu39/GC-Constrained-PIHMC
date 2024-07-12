//
// Created by Jin Bin on 2021/10/02.
//

#ifndef MD_NVE_LJ_CELL_H
#define MD_NVE_LJ_CELL_H

#include "box.h"
#include "matrix3.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
using namespace std;

class cell
{
public:
    string fn;
    int ntype = Consts.ntype;
    // atom *Atoms; // save coordinate after scaled
    Matrix3 lattice_vector; // save lattice vector after scaled
    element *elements;
    double lattice;
    int N_atoms;
    string global_out_dir;
    string *atom_label;
    double *atom_mass;
    string Coordinate;
    double omega = 0.0;
    double electron_number = 0;
    explicit cell(const string &atom_file="", const string &global_out_dir="./");
    ~cell();
    void setup_cell();
    void print_cell_xyz(const string &fn) const;
    void print_stru_file(const string &fn, const int &type=1) const;
private:
    // atom *atoms;
    bool gce = false;
    string *pseudo_fn;
    string *numer_orb;
    bool lcao = false;
    double lat0 = 0.0;
    double lat0_angstrom = 0.0;
    Matrix3 latvec;
    Vector3<double> a1, a2, a3;
    // void scale_forward();
    // void scale_backward();
    void read_atom_species(ifstream &ifa);
    bool read_atom_positions(ifstream &ifpos);
    void scale_lattice();
    static bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, bool restart=1); //NOLINT
    template <class T>
    static void READ_VALUE(ifstream &ifs, T &v)
    {
        ifs >> v;
        ifs.ignore(150, '\n');
    }
    static inline void Cartesian_to_Direct
            (
                    const double &cx,const double &cy,const double &cz,
                    const double &R11,const double &R12,const double &R13,
                    const double &R21,const double &R22,const double &R23,
                    const double &R31,const double &R32,const double &R33,
                    double &dx,double &dy,double &dz)
    {
        static Matrix3 lattice_vector, inv_lat;
        lattice_vector.e11 = R11;
        lattice_vector.e12 = R12;
        lattice_vector.e13 = R13;
        lattice_vector.e21 = R21;
        lattice_vector.e22 = R22;
        lattice_vector.e23 = R23;
        lattice_vector.e31 = R31;
        lattice_vector.e32 = R32;
        lattice_vector.e33 = R33;

        inv_lat = lattice_vector.Inverse();

        static Vector3<double> direct_vec, cartesian_vec;
        cartesian_vec.x = cx;
        cartesian_vec.y = cy;
        cartesian_vec.z = cz;
        direct_vec = cartesian_vec * inv_lat;
        dx = direct_vec.x;
        dy = direct_vec.y;
        dz = direct_vec.z;
    }
};


#endif //MD_NVE_LJ_CELL_H
