# __CPIHMC__
Documents for __CPIHMC__ code version _0.1.0alpha_.

# __Functions__
MD, HMC, CHMC simulations with different types of potential under a _canonical_ or _grand canonical_ ensemble condition. Besides, nuclear quantum effects can be considered through the path integral formalism involved in the method.

1. Potential energy and forces can be calculated by _DP_, _VASP_ and _ABACUS_. (Only _DP_ is supported for the moment.)
2. For the _grand canonical_ ensemble, only CHMC, CPIHMC is supported.
3. Model deviation is only available for the potential energy and forces calculated by _DP_. Meanwhile, all situations can be performed in a model deviation mode.

# __Units__
### __Physical Quantities__
- __mass:__ u, $1$ u = $1.66053886\times 10^{-27}$ kg

- __length:__ Bohr, $1$ Bohr = $0.529177\AA$

- __time:__ fs

- __energy:__ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$, $1$ J = $2.15055\times 10^{17}$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$, $1$ eV = $0.0344556$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$

- __force:__ u $\cdot$ Bohr $\cdot$ fs $^{-2}$, $1$ eV/$\AA$ = $0.0182331$ u $\cdot$ Bohr $\cdot$ fs $^{-2}$
### __Physical Constants__
- $k$ = $2.96915\times 10^{-6}$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}\cdot$ K $^{-1}$

- $h$ = $0.142497$ u $\cdot$ Bohr $^2\cdot$ fs $^{-1}$

- $\hbar$ = $0.0226791$ u $\cdot$ Bohr $^2\cdot$ fs $^{-1}\cdot$ rad $^{-1}$

# __Input Files__
## `INPUT`
Simulation parameters should be set in __`INPUT`__, and the parameters not in the __`INPUT`__ are set to be default values. All optional parameters are listed later in this document.

## `STRU` (, `BEADS`)
__`STRU`__ contains the initial structure in the format of _ABACUS_. However, if necessary, the electron number of the system should be put in the beginning of __`STRU`__ by a title __ELECTRON_NUMBER__ with a following line of the electron number value.

As for __path integral__, if there's a need to set different coordinates of beads at the beginning, set the centroid coordinates in __`STRU`__, and put different coordinates of beads in __`BEADS`__. For this situation, don't forget to add the parameter __Beads_File__ into __`INPUT`__.

> __Attention__: The unit of coordinates in __`BEADS`__ is Bohr.

# __(DP potential)__
For the _DP_ situation, don't forget to prepare a DP potential and refer to the potential with the parameter __Deep_Pot_Model__ in __`INPUT`__.

# Output Files
### __`ALL_INPUT`__
The value of all parameters used in the simulation.

### __`PHY_QUANT`__
Selected physical quantities for each several step(s). The output interval can be set by the paramter __Phy_Quant_Intvl__ in __`INPUT`__.

### __`ALL_STRU`__
Dumped structures along the simulation process.

### __(`MODEL_DEVI`)__
Model deviation results of virials and forces between different _DP_ potentials. Only output in the __model deviation__ mode.

# Parameters
### __General Settings__
- __Stru_File__
default: STRU (optional)

> the structure file. It's necessary only when you want to rename the structure file.

- __Pot_Type__
default: DP (optional)

> the style used to obtain the potential energy and forces.
> 
> choices: DP

- __Deep_Pot_Model__
default: (optional)

> the DP potential for the system. It's necessary when the DP potential is used for the potential energy and forces.

- __Simu_Type__
default: CHMC (optional)

> the simulation type.
> 
> choices: MD, HMC, CHMC

- __Steps__
(required)

> the number of simulation steps. A step means a sampling process between two adjacent acceptance judges.

### __System Settings__
- __Temp__
default: 300

> unit: K
>                  
> the temperature of the simulation.

- __N_Type__
(required)

> the number of element(s).

- __Element_Type__
(required)

> different element type(s) at the order in DP potential. The parameters are always required but only used for the order when the parameter __Element_Index__ is not set.

- __Element_Index__
default: / (optional)

> different element index(s) used in DP potential mapping to element type(s) one-to-one. The parameter only displays when it has been set up.

- __Wall__
default: None (optional)

> unit (wall position): Bohr
>
> the hard wall(s) set to limit atom coordinates in the cell. The wall can be set on 6 directions. However, the function is not available in the MD simulations. The format of this parameter is __string-floating point pair(s)__ constituted by a wall name and a number defining the position of the wall, such as _xlo 129.3_.
>
> choices (wall name): xlo, xhi, ylo, yhi, zlo, zhi

- __Group__
default: None (optional)

> the group(s) of atom(s) defined with range(s). The format of this parameter is a group name followed by a range of the index(es) of atom(s) composing the group. The group name can be set up arbitrarily like the name of a variable in the code. The range is defined in the format b:e[:s], which means the atom(s) with the index(es) in the interval [b, e) under a stride s is(are) included in the group. The default stride is 1.
>
> example:
> Group bottom 1:6 Proton 36:40:2
> (In this example, two groups named "bottom" and "Proton" are defined. The group named "bottom" contains five atoms with the indexes 1, 2, 3, 4 and 5. And the group named "Proton" contains two atoms with the indexes 36 and 38.)

### __Molecular Dynamics Settings__
- __Time_Step__
default: 1.0

> unit: fs
>
> the time step of MD simulations.

### __Hybird Monte Carlo Settings__
- __N_Evol_Step__
default: 3 (optional)

> the number of V-V steps between two adjacent judges.

- __Time_Step__
default: 1.0

> unit: fs
> 
> the time step of V-V step in HMC simulations.

- __Mass_Scal__
default: 1.0 (optional)

> the mass scaling coefficient used in the HMC scheme. Decreasing the masses of all elements owns the similar effect of increasing the velocity sampling temperature.

### __Constraint Settings__
- __Hybrid_Monte_Carlo_Ratio__
default: 0.8 (optional)

> the ratio of the HMC part in the centroid move. The centroid move contains HMC part and MC part.

- __Virt Atom__
default: None (optional)

> the virtual atom(s) used in the definition of the reaction coordinate. The virtual atom is the centroid of the atoms composing it.
> The format of the parameter is a virtual atom name followed by the indexes of atoms composing the virtual atom. The virtual atom name can be set up arbitrarily like the name of a variable in the code.
>
> example:
> Virt_Atom Cent_0 9 6 16 Cent_1 1 36
> (In this example, two virtual atoms named "Cent_0" and "Cent_1" are defined. The virtual atom named "Cent_0" is the centroid of three atoms with the indexes 9, 6 and 16. And the virtual atom named "Cent_1" is the centroid of two atoms with the indexes 1 and 36.)

- __Rxn_Coord__
default: None (optional)

> unit (reaction coordinate value): Bohr
> 
> the reaction coordinate(s) used in the constrained simulations. The format of the parameter is a reaction coordinate type followed by three types of parameters for this reaction coordinate, the indexes of atoms involved in the reaction coordinate, the reaction coordinate value and the parameters for the MC part to evolve the reaction coordinate. The atom index can be either the index of the real atom or the name of the virtual atom.
> 
> choices (reaction coordinate type): DIST, DIFF
>
> __DIST__: the distance between two particles. Two atom indexes are required in the order $i_1$ $i_2$, then the reaction coordinate is $|\vec{r}_{i_1}-\vec{r}_{i_2}|$. Only one parameter is required for the MC part, and the parameter denotes the radius of the random sphere.
> 
> example:
> DIST Cent_0 1 0.6 0.1
> (In this example, the reaction coordinate is the distance between the virtual atom named "Cent_0" and the real atom with the index 1. The reaction coordinate value is 0.6 Bohr, and the radius of the random sphere used in the MC part is 0.1.)
>
> __DIFF__: the difference of particles' distances. Three atom indexes are required in the order $i_1$ $i_2$ $i_3$, then the reaction coordinate is $|\vec{r}_{i_2}-\vec{r}_{i_1}|-|\vec{r}_{i_2}-\vec{r}_{i_3}|$. Three parameters are required for the MC part, the random widths of the length, the theta and the phi.
>
> example:
> DIFF 6 Cent_1 9 -0.1 0.4 0.1309 1.0472
> (In this example, the reaction coordinate is the difference of the distances between the virtual atom named "Cent_1" with two real atoms with indexes 6 and 9 respectively. The reaction coordinate value is -0.1 Bohr, and three parameters used in the MC part are 0.4, 0.1309 and 1.0472.)

### __Grand Canonical Ensemble Settings__
- __Elec_Num_Ratio__
default: 0.0 (optional)

> the ratio of the electron number variation part in all moves. The electron number variation part will only be involved in the simulation if this parameter doesn't equal to 0.0.

- __Mu__
default: 0.0 (optional)

> unit: eV
> 
> the electrochemical potential for the _grand canonical_ ensemble.

- __Elec_Num_Range__
default: 0.0 0.0 (optional)

> unit: 1
>
> the range of the electron number can be sampled. Two numbers of the parameter are the lower and higher bounds of the range, respectively. The electron number of the initial structure should locate in this range.

- __Elec_Num_Width__
default: 0.04 (optional)

> unit: 1
> 
> the change width of the electron number.

### __Path Integral Settings__
- __Beads_File__
default: (optional)

> the initial beads coordinates file. If the parameter isn't given, all beads are the same as the centroid at the beginning.

- __N_Bead__
default: 1 (optional)

> the number of beads.

- __N_Change_Bead__
default: 1 (optional)

> the number of bead(s) changed in each step. This parameter should be smaller than the parameter __N_Bead__.

- __Bead_Index__
default: None (optional)

> the index(s) of different bead atom(s) in coordinates.The atom index can be either the index of the real atom or the name of the group.

### __Model Deviation Settings__
- __Model_Devi_Deep_Pot_Models__
default: None (optional)

> the DP potentials used for the model deviations. More than one DP potentials should be provided, and the DP potential used for the simulation isn't contained inherently.

- __Model_Devi_Intvl__
default: 100 (optional)
> the interval between two outputs of model deviation results.

- __Model_Devi_File__
default: MODEL_DEVI (optional)

> the file to output the model deviation results.

### __Output Settings__
- __Phy_Quant_File__
default: PHY_QUANT (optional)
                    
> the file to output different physical quantities.

- __Col_Width__
default: 12 (optional)

> the width of the column in the file to output different physical quantities set by the parameter __Phy_Quant_File__.

- __Out_Phy_Quant__
default: ke pe etotal (optional)

> the physical quantities to be output in the file set by __Phy_Quant_File__. The physical quantities can be chosen contain kinetic energy, potential energy, total energy, internal energy, electron number, reaction coordinate value, mean force on the left atom, mean force on the right atom, mean force, accumulative average of potential energy.
>
> choices: ke pe etotal eint ne rc mfl mfr mf pe_ave
> (if more than one reaction coordinates are defined in the simulation, the reaction coordinate value, mean force of each reaction coordinate are labeled by a suffix index beginning with 0. For example, the reaction coordinate value of the second reaction coordinate is _rc_1_)

- __Phy_Quant_Digits__
default: ke 6 pe 6 etotal 6 (optional)

> the decimal digits of physical quantity(ies) to be output in the file set by __Phy_Quant_File__. The format of this parameter is string-integer pair(s) constituted by a physical quantity and a number defining the decimal digit of the physical quantity, such as _ke 2_.

- __Phy_Quant_Intvl__
default: 1 (optional)

> the number of steps in the interval to output the physical quantity(ies).

- __Stru_Intvl__
default: 100 (optional)

> the number of steps in the interval between two outputs of the structures.