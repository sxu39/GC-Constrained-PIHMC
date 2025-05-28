#include "cell.h"
using namespace cpihmc;
using namespace std;

cpihmc::cell::cell(const input &Input):AtomLabel(nullptr), AtomMass(nullptr), PseudoFileName(nullptr), OrbitalFileName(nullptr), 
                                       Magnetization(nullptr), AtomMagnet(nullptr), AtomMagnetVec(nullptr), SetElecNum(false), ElecNum(0.0), 
                                       LatConst(0.0), LatConstAng(0.0), Omega(0.0), NType(Input.get_n_type()), NAtoms(0)
{
    const string &StruFile = Input.get_stru_file();
    if (StruFile.size()) setup_cell(StruFile);
}

cpihmc::cell::~cell()
{
    if (AtomLabel) delete [] AtomLabel;
    if (AtomMass) delete [] AtomMass;
    if (PseudoFileName) delete [] PseudoFileName;
    if (OrbitalFileName) delete [] OrbitalFileName;
    if (Magnetization) delete [] Magnetization;
    if (AtomMagnet)
    {
        for (index_t i = 0 ; i < NType ; ++i)
            delete [] AtomMagnet[i];
        delete [] AtomMagnet;
    }
    if (AtomMagnetVec)
    {
        for (index_t i = 0 ; i < NType ; ++i)
            delete [] AtomMagnetVec[i];
        delete [] AtomMagnetVec;
    }
}

void cpihmc::cell::setup_cell(const string &StruFile)
{
	assert(NType>0);

	bool_t ok = true;
	bool_t ok2 = true;

    // open structure file.
    ifstream InFileStru(StruFile, ios::in);
    if (!InFileStru) ok = false;

    if(ok)
    {
        read_atom_species(InFileStru);
        ok2 = read_atom_positions(InFileStru);
    }
	
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, LatVec must be read in.
	//========================================================
	assert(LatConst > 0.0);
	Omega = abs(LatVec.det()) * pow(LatConst, 3);
}

void cpihmc::cell::read_atom_species(ifstream &InFileStru)
{
    if (AtomLabel) delete [] AtomLabel;
    if (AtomMass) delete [] AtomMass;
    if (PseudoFileName) delete [] PseudoFileName;
    AtomLabel = new string [NType]; // atom labels
    AtomMass = new prec_t [NType]; //atom masses
    PseudoFileName = new string [NType]; //file name of pseudopotential

    //==========================================
    // read in the electron number if necessary
    //==========================================
    if(scan_begin(InFileStru, "ELECTRON_NUMBER"))
    {
        read_value(InFileStru, ElecNum);
        SetElecNum = true;
    }

    //==========================================
    // read in information of each type of atom
    //==========================================
    if(scan_begin(InFileStru, "ATOMIC_SPECIES"))
    {    
        InFileStru.ignore(300, '\n');
        for (index_t i = 0 ; i < NType ; ++i)
        {
            string OneLine, OneStr;
            getline(InFileStru, OneLine);
            stringstream StrStr;
            StrStr << OneLine;
            StrStr >> AtomLabel[i] >> AtomMass[i];
            PseudoFileName[i] = "auto";
            if (StrStr >> OneStr)
                if (OneStr[0] != '#')
                    PseudoFileName[i] = OneStr;
        }
    }

    //==============================================
    // read in numerical orbital files if necessary
    //==============================================
    if(scan_begin(InFileStru, "NUMERICAL_ORBITAL"))
    {
        if (OrbitalFileName) delete [] OrbitalFileName;
        OrbitalFileName = new string [NType];
        for (index_t i = 0 ; i < NType ; ++i)
            InFileStru >> OrbitalFileName[i];
    }

    //==========================
    // read in lattice constant
    //==========================
    if(scan_begin(InFileStru, "LATTICE_CONSTANT"))
    {
        read_value(InFileStru, LatConst);
        LatConstAng = LatConst * Bohr2Ang;
    }

    //===========================
    // Read in lattice vectors
    //===========================
    if(scan_begin(InFileStru, "LATTICE_VECTORS"))
    {
        // Reading lattice vectors. notice here that only one cpu read these parameters.
        InFileStru >> LatVec.e11 >> LatVec.e12;
        read_value(InFileStru, LatVec.e13);
        InFileStru >> LatVec.e21 >> LatVec.e22;
        read_value(InFileStru, LatVec.e23);
        InFileStru >> LatVec.e31 >> LatVec.e32;
        read_value(InFileStru, LatVec.e33);
    }
} // end read_atom_species

const bool_t cpihmc::cell::read_atom_positions(ifstream &InFileStru)
{
    stringstream StrWarn;
    if(scan_begin(InFileStru, "ATOMIC_POSITIONS"))
    {
        read_value(InFileStru, CoordType);
        if (CoordType != "Cartesian" 
            && CoordType != "Direct" 
            && CoordType != "Cartesian_angstrom"
            && CoordType != "Cartesian_au"
            && CoordType != "Cartesian_angstrom_center_xy"
            && CoordType != "Cartesian_angstrom_center_xz"
            && CoordType != "Cartesian_angstrom_center_yz"
            && CoordType != "Cartesian_angstrom_center_xyz"
            )
        {
            StrWarn << " There are several options for you:" << endl;
            StrWarn << " Direct" << endl;
            StrWarn << " Cartesian_angstrom" << endl;
            StrWarn << " Cartesian_au" << endl;
            StrWarn << " Cartesian_angstrom_center_xy" << endl;
            StrWarn << " Cartesian_angstrom_center_xz" << endl;
            StrWarn << " Cartesian_angstrom_center_yz" << endl;
            StrWarn << " Cartesian_angstrom_center_xyz" << endl;
            cerr << StrWarn.str();
            return false; // means something wrong
        }

        vec3<prec_t> Vec;
        vec3<bool_t> MoveVec;
        size_t NAtom = 0;
        NAtoms = 0;

        //======================================
        // calculate total number of atoms
        // and adjust the order of atom species
        //======================================
        assert(NType>0);
        if (Magnetization) delete [] Magnetization;
        Magnetization = new prec_t [NType];
        if (AtomMagnet)
        {
            for (index_t i = 0 ; i < NType ; ++i)
                delete [] AtomMagnet[i];
            delete [] AtomMagnet;
        }
        AtomMagnet = new prec_t *[NType];
        if (AtomMagnetVec)
        {
            for (index_t i = 0 ; i < NType ; ++i)
                delete [] AtomMagnetVec[i];
            delete [] AtomMagnetVec;
        }
        AtomMagnetVec = new vec3<prec_t> *[NType];
        string Label;
        index_t AtomIndex = 0;
        for (index_t it = 0 ; it < NType ; ++it)
        {   
            //=======================================
            // read in atom label
            // start magnetization
            //=======================================
            read_value(InFileStru, Label);
            if (Label == AtomLabel[it]) Elements.push_back(element{Label, AtomMass[it], Atoms});
            else
            {
                StrWarn << " Label orders in ATOMIC_POSITIONS and ATOMIC_SPECIES sections do not match!" << endl;
                StrWarn << " Label read from ATOMIC_POSITIONS is " << Elements[it].Label << endl;
                StrWarn << " Label from ATOMIC_SPECIES is " << AtomLabel[it] << endl;
                cerr << StrWarn.str();
                return false;
            }

            read_value(InFileStru, Magnetization[it]);

            //=========================
            // read in atom number
            //=========================
            read_value(InFileStru, NAtom);
            Elements[it].Number = NAtom;
            NAtoms += NAtom;
            Atoms.resize(NAtoms, atom{Elements[it].Mass});
            Elements[it].set_begin(NAtoms-NAtom);

            if (NAtom < 0) return false;
            if (NAtom > 0)
            {
                AtomMagnet[it] = new prec_t [NAtom];
                AtomMagnetVec[it] = new vec3<prec_t> [NAtom];
                for (index_t ia = 0 ; ia < NAtom ; ++ia)
                {
                    InFileStru >> Vec.x >> Vec.y >> Vec.z;
                    MoveVec.x = true;
                    MoveVec.y = true;
                    MoveVec.z = true ;

                    string Tmpid;
                    Tmpid = InFileStru.get();

                    if((index_t)Tmpid[0] < 0)
                    {
                        StrWarn << "read_atom_positions, mismatch in atom number for atom type: " << Elements[it].Label << endl;
                        cerr << StrWarn.str();
                        exit(1); 
                    }

                    // read if catch goodbit before "\n" and "#"
                    while ((Tmpid != "\n") && (InFileStru.good()) && (Tmpid !="#") )
                    {
                        Tmpid = InFileStru.get();
                        // old method of reading frozen ions
                        char Tmp = (char)Tmpid[0];
                        if (Tmp >= 48 && Tmp <= 57)
                        {
                            MoveVec.x = stoi(Tmpid);
                            InFileStru >> MoveVec.y >> MoveVec.z ;
                        }
                        // new method of reading frozen ions and velocities
                        if (Tmp >= 'a' && Tmp <= 'z')
                        {
                            InFileStru.putback(Tmp);
                            InFileStru >> Tmpid;
                        }
                        if (Tmpid == "m")
                            InFileStru >> MoveVec.x >> MoveVec.y >> MoveVec.z;
                        else if (Tmpid == "v" || Tmpid == "vel" || Tmpid == "velocity")
                            InFileStru >> Atoms[AtomIndex].Vel.x >> Atoms[AtomIndex].Vel.y >> Atoms[AtomIndex].Vel.z;
                        else if (Tmpid == "mask") Masks.push_back(AtomIndex);
                        else if (Tmpid == "mag" || Tmpid == "magmom")
                        {
                            prec_t Tmpamg = 0;
                            InFileStru >> Tmpamg;
                            Tmp = InFileStru.get();
                            while (Tmp==' ')
                                Tmp=InFileStru.get();
                            
                            if ((Tmp >= 48 && Tmp <= 57) || Tmp=='-')
                            {
                                InFileStru.putback(Tmp);
                                InFileStru >> AtomMagnetVec[it][ia].y >> AtomMagnetVec[it][ia].z;
                                AtomMagnetVec[it][ia].x = Tmpamg;
                            }
                            else
                            {
                                InFileStru.putback(Tmp);
                                AtomMagnet[it][ia] = Tmpamg;
                            }
                        }
                        else if (Tmpid == "angle1")
                        {
                            prec_t Angle1;
                            InFileStru >> Angle1;
                        }
                        else if (Tmpid == "angle2")
                        {
                            prec_t Angle2;
                            InFileStru >> Angle2;
                        }    
                    }
                    // move to next line
                    while ((Tmpid != "\n") && (InFileStru.good()))
                        Tmpid = InFileStru.get();
            
                    if (CoordType == "Direct")
                    {
                        // change Vec from Direct to Cartesian
                        Atoms[AtomIndex].Coord = Vec * LatVec;
                    }
                    else if (CoordType == "Cartesian")
                        Atoms[AtomIndex].Coord = Vec; // in unit LatConst
                    else if (CoordType == "Cartesian_angstrom")
                        Atoms[AtomIndex].Coord = Vec / LatConstAng;
                    else if (CoordType == "Cartesian_angstrom_center_xy")
                    {
                        // calculate lattice center
                        vec3<prec_t> LatCenter{(LatVec.e11+LatVec.e21+LatVec.e31)/2.0, (LatVec.e12+LatVec.e22+LatVec.e32)/2.0, 0.0};
                        Atoms[AtomIndex].Coord = Vec / LatConstAng + LatCenter;
                    }
                    else if (CoordType == "Cartesian_angstrom_center_xz")
                    {
                        // calculate lattice center
                        vec3<prec_t> LatCenter{(LatVec.e11+LatVec.e21+LatVec.e31)/2.0, 0.0, (LatVec.e13+LatVec.e23+LatVec.e33)/2.0};
                        Atoms[AtomIndex].Coord = Vec / LatConstAng + LatCenter;
                    }
                    else if (CoordType == "Cartesian_angstrom_center_yz")
                    {
                        // calculate lattice center
                        vec3<prec_t> LatCenter{0.0, (LatVec.e12+LatVec.e22+LatVec.e32)/2.0, (LatVec.e13+LatVec.e23+LatVec.e33)/2.0};
                        Atoms[AtomIndex].Coord = Vec / LatConstAng + LatCenter;
                    }
                    else if (CoordType == "Cartesian_angstrom_center_xyz")
                    {
                        // calculate lattice center
                        vec3<prec_t> LatCenter{(LatVec.e11+LatVec.e21+LatVec.e31)/2.0, (LatVec.e12+LatVec.e22+LatVec.e32)/2.0, 
                                               (LatVec.e13+LatVec.e23+LatVec.e33)/2.0};
                        Atoms[AtomIndex].Coord = Vec / LatConstAng + LatCenter;
                    }
                    else if (CoordType == "Cartesian_au")
                        Atoms[AtomIndex].Coord = Vec / LatConst;
                    
                    Atoms[AtomIndex].Move = MoveVec;
                    ++AtomIndex;
                } // end for NAtom
            } // end NAtom
        } // end for NType
    } // end scan_begin

    for (index_t it = 0 ; it < NType ; ++it)
        for (index_t ia = 0 ; ia < Elements[it].Number ; ++ia)
            Elements[it][ia].Coord *= LatConst;
    return true;
} // end read_atom_positions