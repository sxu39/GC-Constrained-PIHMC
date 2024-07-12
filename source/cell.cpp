//
// Created by Jin Bin on 2022/02/06.
//

#include "cell.h"

cell::cell(const string &atom_file, const string &global_out_dir):fn(atom_file), global_out_dir(global_out_dir){
    if (strcmp(atom_file.c_str(), ""))
        setup_cell();
}

cell::~cell(){
    delete [] elements;
    delete [] atom_label;
    delete [] atom_mass;
    delete [] pseudo_fn;
    delete [] numer_orb;
}

void cell::setup_cell(){
    assert(ntype>0);
    elements = new element [ntype];

    bool ok = true;
    bool ok2 = true;

    ifstream ifa(fn.c_str(), ios::in);
    if (!ifa)
        ok = false;
    if(ok)
    {
        this->read_atom_species(ifa);

        //==========================
        // call read_atom_positions
        //==========================
        ok2 = this->read_atom_positions(ifa);
    }
    //========================================================
    // Calculate unit cell volume
    // the reason to calculate volume here is
    // Firstly, latvec must be read in.
    //========================================================
    assert(lat0 > 0.0);
    this->omega = abs( latvec.Det()) * this->lat0 * lat0 * lat0;
}

void cell::read_atom_species(ifstream &ifa){
    pseudo_fn = new string [ntype];
    atom_label = new string [ntype];
    atom_mass = new double [ntype];
    //==========================================
    // read in the electron number if necessary
    //==========================================
    if( SCAN_BEGIN(ifa, "ELECTRON_NUMBER") ){
        READ_VALUE(ifa, electron_number);
        gce = true;
    }

    //==========================================
    // read in information of each type of atom
    //==========================================
    if( SCAN_BEGIN(ifa, "ATOMIC_SPECIES") ){
        for (int i = 0 ; i < ntype ; ++i){
            ifa >> atom_label[i] >> atom_mass[i];
            READ_VALUE(ifa, pseudo_fn[i]);
        }
    }

    //==============================================
    // read in numerical orbital files if necessary
    //==============================================
    if (SCAN_BEGIN(ifa, "NUMERICAL_ORBITAL")){
        lcao = true;
        numer_orb = new string [ntype];
        for (int i = 0 ; i < ntype ; ++i)
            READ_VALUE(ifa, numer_orb[i]);
    }

    //==========================
    // read in lattice constant
    //==========================
    if( SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
    {
        READ_VALUE(ifa, lat0);
        lat0_angstrom = lat0 * 0.529177;
    }

    //===========================
    // Read in latticies vector
    //===========================
    if( SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
    {
        // Reading lattice vectors. notice
        // here that only one cpu read these
        // parameters.
        ifa >> latvec.e11 >> latvec.e12;
        READ_VALUE(ifa, latvec.e13);
        ifa >> latvec.e21 >> latvec.e22;
        READ_VALUE(ifa, latvec.e23);
        ifa >> latvec.e31 >> latvec.e32;
        READ_VALUE(ifa, latvec.e33);

        // lattice vectors in another form.
        a1.x = latvec.e11;
        a1.y = latvec.e12;
        a1.z = latvec.e13;

        a2.x = latvec.e21;
        a2.y = latvec.e22;
        a2.z = latvec.e23;

        a3.x = latvec.e31;
        a3.y = latvec.e32;
        a3.z = latvec.e33;

    }
    lattice = lat0 * latvec.e11;
    this->scale_lattice();
}

bool cell::read_atom_positions(ifstream &ifpos)
{
    bool use_xyz = false;

    if( SCAN_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        READ_VALUE(ifpos, Coordinate);
        if(Coordinate != "Cartesian"
           && Coordinate != "Direct"
           && Coordinate != "Cartesian_angstrom"
           && Coordinate != "Cartesian_au"
           && Coordinate != "Cartesian_angstrom_center_xy"
           && Coordinate != "Cartesian_angstrom_center_xz"
           && Coordinate != "Cartesian_angstrom_center_yz"
           && Coordinate != "Cartesian_angstrom_center_xyz"
                )
            return 0;

        Vector3<double> v;
        Vector3<int> mv;
        int na = 0;
        this->N_atoms = 0;

        //======================================
        // calculate total number of atoms
        // and adjust the order of atom species
        //======================================
        assert(ntype>0);
        for (int it = 0 ; it < ntype ; ++it){
            string mag;
            //=======================================
            // (1) read in atom label
            // start magnetization
            //=======================================
            READ_VALUE(ifpos, elements[it].label);
            bool found = false;
            for (int it2 = 0 ; it2 < ntype ; ++it2){
                if(elements[it].label == this->atom_label[it]){
                    found = true;
                    // scale the mass for HMC evolution, default mass scaling coefficient is 1
                    elements[it].mass = this->atom_mass[it] * Consts.M_scaling;
                    break;
                }
            }
            if(!found)
                return 0;
            READ_VALUE(ifpos, mag);

            //=========================
            // (3) read in atom number
            //=========================
            READ_VALUE(ifpos, na);

            this->N_atoms += na;
            if (na <= 0)
                return 0;
            if (na > 0)
            {
                this->elements[it].number = na;
                this->elements[it].atoms = new atom [N_atoms];
                atom * &atoms = this->elements[it].atoms;
                for (int ia = 0;ia < na; ia++)
                {
                    if(use_xyz)
                    {
                        string tmpid;
                        ifpos >> tmpid;
                        if(tmpid != elements[it].label)
                            return 0;
                        else
                        {
                            ifpos >> v.x >> v.y >> v.z;

                            mv.x = true;
                            mv.y = true;
                            mv.z = true;
                        }
                    }
                    else
                    {
                        ifpos >> v.x >> v.y >> v.z
                              >> mv.x >> mv.y >> mv.z;
                    }

                    ifpos.ignore(150, '\n');

                    atoms[ia].move.x = mv.x;
                    atoms[ia].move.y = mv.y;
                    atoms[ia].move.z = mv.z;

                    if(Coordinate=="Direct")
                    {
                        // change v from direct to cartesian,
                        // the unit is pw.lat0
                        Vector3<double> v_c = v * latvec;
                        atoms[ia].r.x = v_c.x;
                        atoms[ia].r.y = v_c.y;
                        atoms[ia].r.z = v_c.z;
                    }
                    else if(Coordinate=="Cartesian")
                    {
                        atoms[ia].r.x = v.x;
                        atoms[ia].r.y = v.y;
                        atoms[ia].r.z = v.z;
                    }
                    else if(Coordinate=="Cartesian_angstrom")
                    {
                        v = v / 0.529177 / lat0;
                        atoms[ia].r.x = v.x;
                        atoms[ia].r.y = v.y;
                        atoms[ia].r.z = v.z;
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xy")
                    {
                        // calculate lattice center
                        v = v / 0.529177 / lat0;
                        atoms[ia].r.x = v.x + (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        atoms[ia].r.y = v.y + (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        atoms[ia].r.z = v.z;
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xz")
                    {
                        // calculate lattice center
                        v = v / 0.529177 / lat0;
                        atoms[ia].r.x = v.x + (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        atoms[ia].r.y = v.y;
                        atoms[ia].r.z = v.z + (latvec.e13 + latvec.e23 + latvec.e33)/2.0;
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_yz")
                    {
                        // calculate lattice center
                        v = v / 0.529177 / lat0;
                        atoms[ia].r.x = v.x;
                        atoms[ia].r.y = v.y + (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        atoms[ia].r.z = v.z + (latvec.e13 + latvec.e23 + latvec.e33)/2.0;
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xyz")
                    {
                        // calculate lattice center
                        v = v / 0.529177 / lat0;
                        atoms[ia].r.x = v.x + (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        atoms[ia].r.y = v.y + (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        atoms[ia].r.z = v.z + (latvec.e13 + latvec.e23 + latvec.e33)/2.0;
                    }
                    else if(Coordinate=="Cartesian_au")
                    {
                        v = v / lat0;
                        atoms[ia].r.x = v.x;
                        atoms[ia].r.y = v.y;
                        atoms[ia].r.z = v.z;
                    }
                }
            }
        }

    }
    // this->scale_forward();
    for (int i = 0 ; i < ntype ; ++i){
        for (int j = 0 ; j < this->elements[i].number ; ++j){
            this->elements[i].atoms[j].r.x *= lat0;
            this->elements[i].atoms[j].r.y *= lat0;
            this->elements[i].atoms[j].r.z *= lat0;
        }
    }
    // this->print_cell_xyz("STRU_READIN_ADJUST.xyz");
    return 1;
}//end read_atom_positions

void cell::print_cell_xyz(const string &fn) const
{
    stringstream ss;
    ss << global_out_dir << fn;

    ofstream ofs(ss.str().c_str());

    ofs << N_atoms << endl;
    ofs << "None" << endl;
    ofs << atom_label << endl;
    for (int it = 0 ; it < ntype ; ++it){
        for (int ia=0; ia < elements[it].number ; ia++)
        {
            ofs << atom_label
                << " " << elements[it].atoms[ia].r.x * 0.529177
                << " " << elements[it].atoms[ia].r.y * 0.529177
                << " " << elements[it].atoms[ia].r.z * 0.529177 << endl;
        }
    }

    ofs.close();
}

void cell::print_stru_file(const string &fn, const int &type) const
{
    ofstream ofs(fn.c_str());

    if (gce){
        ofs << "ELECTRON_NUMBER\n" << electron_number << endl;
        ofs << endl;
    }

    ofs << "ATOMIC_SPECIES" << endl;
    ofs << setprecision(12);

    for(int it=0; it<ntype; it++)
    {
        //modified by zhengdy 2015-07-24
        ofs << atom_label[it] << " " << atom_mass[it] << " " << pseudo_fn[it] << endl;
    }

    if (lcao){
        ofs << "\nNUMERICAL_ORBITAL" << endl;
        for (int it = 0 ; it < ntype ; it++)
            ofs << numer_orb[it] << endl;
    }

    ofs << "\nLATTICE_CONSTANT" << endl;
    //modified by zhengdy 2015-07-24
    ofs << lat0 << endl;

    ofs << "\nLATTICE_VECTORS" << endl;
    ofs << latvec.e11 << " " << latvec.e12 << " " << latvec.e13 << " #latvec1" << endl;
    ofs << latvec.e21 << " " << latvec.e22 << " " << latvec.e23 << " #latvec2" << endl;
    ofs << latvec.e31 << " " << latvec.e32 << " " << latvec.e33 << " #latvec3" << endl;

    ofs << "\nATOMIC_POSITIONS" << endl;

    if(type==1)
    {
        ofs << "Cartesian" << endl;
        for(int it=0; it<ntype; it++)
        {
            ofs << endl;
            ofs << elements[it].label << " #label" << endl;
            ofs << 0 << " #magnetism" << endl;
            //2015-05-07, modify
            //ofs << atoms[it].nwl << " #max angular momentum" << endl;
            //xiaohui modify 2015-03-15
            //for(int l=0; l<=atoms[it].nwl; l++)
            //{
            //	ofs << atoms[it].l_nchi[l] << " #number of zeta for l=" << l << endl;
            //}
            ofs << elements[it].number << " #number of atoms" << endl;
            for(int ia=0; ia<elements[it].number; ia++)
            {
                ofs << elements[it].atoms[ia].r.x / lat0 << " "
                    << elements[it].atoms[ia].r.y / lat0 << " "
                    << elements[it].atoms[ia].r.z / lat0 << " "
                    << elements[it].atoms[ia].move.x << " "
                    << elements[it].atoms[ia].move.y << " "
                    << elements[it].atoms[ia].move.z << endl;
            }
        }
    }
    else if(type==2)
    {
        ofs << "Direct" << endl;
        for (int it = 0 ; it < ntype ; ++it){
            ofs << endl;
            ofs << elements[it].label << " #label" << endl;
            ofs << 0 << " #magnetism" << endl;
            //ofs << atoms[it].nwl << " #max angular momentum" << endl;
            //xiaohui modify 2015-03-15
            //for(int l=0; l<=atoms[it].nwl; l++)
            //{
            //	ofs << atoms[it].l_nchi[l] << " #number of zeta for l=" << l << endl;
            //}
            ofs << elements[it].number << " #number of atoms" << endl;
            for(int ia=0; ia<elements[it].number; ia++)
            {
                double dx,dy,dz;
                Cartesian_to_Direct(elements[it].atoms[ia].r.x / lat0,
                                    elements[it].atoms[ia].r.y / lat0,
                                    elements[it].atoms[ia].r.z / lat0,
                                    latvec.e11, latvec.e12, latvec.e13,
                                    latvec.e21, latvec.e22, latvec.e23,
                                    latvec.e31, latvec.e32, latvec.e33,
                                    dx,dy,dz);
                ofs << dx << " " << dy << " " << dz << " "
                    << elements[it].atoms[ia].move.x << " "
                    << elements[it].atoms[ia].move.y << " "
                    << elements[it].atoms[ia].move.z << endl;
            }
        }
    }

    ofs.close();
}

void cell::scale_lattice(){
    lattice_vector = latvec * lat0;
}

bool cell::SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart)
{
    string SearchName;
    bool find = false;
    if (restart)
    {
        ifs.clear();
        ifs.seekg(0);
    }
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> SearchName;
        if ( SearchName == TargetName)
        {
            find = true;
            break;
        }
    }
    return find;
}