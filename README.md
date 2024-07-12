# __CPIHMC__
Documents for __CPIHMC__ code version _0.1.8beta_.

# __Functions__
MD, HMC, CHMC, CMC simulations with different types of potential under a _canonical_ or _grand canonical_ ensemble condition. Besides, nuclear quantum effects can be considered through the path integral involved in the method.

1. Potential energy and forces can be calculated by _DP_, _VASP_ and _ABACUS_.

2. For "__difference__" reaction coordinate, CHMC, CMC, CPIHMC, CPIMC methods are supported.

3. For "__distance__" reaction coordinate, CHMC, CMC, CPIHMC methods are supported.
4. For the _grand canonical_ ensemble, only HMC, CHMC, CPIHMC is supported.
5. Model deviation is only available for the potential energy and forces calculated by _DP_. Meanwhile, except for the simulation under a _grand canonical_ ensemble condition using the discrete DP potentials, all situations can be performed in a model deviation mode.

# __Units__
### __Physical Quantities__
- __mass:__ u, $1$ u = $1.66053886\times 10^{-27}$ kg

- __length:__ Bohr, $1$ Bohr = $0.529\AA$

- __time:__ fs

- __energy:__ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$, $1$ J = $2.1520\times 10^{17}$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$, $1$ eV = $0.03448$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}$

- __force:__ u $\cdot$ Bohr $\cdot$ fs $^{-2}$, $1$ eV/$\AA$ = $0.01824$ u $\cdot$ Bohr $\cdot$ fs $^{-2}$
### __Physical Constants__
- _k_ = $2.97116\times 10^{-6}$ u $\cdot$ Bohr $^2\cdot$ fs $^{-2}\cdot$ K $^{-1}$

- _h_ = $0.142612$ u $\cdot$ Bohr $^2\cdot$ fs $^{-1}$

# __Input Files__
## `INPUT`
Simulation parameters should be set in __`INPUT`__, and the parameters not in the __`INPUT`__ are set to be default values. All optional parameters are listed later in this document.

## __`STRU` (, `BEADS`)__
__`STRU`__ contains the initial structure in the format of _ABACUS_. However, if necessary, the electron number of the system should be put in the beginning of __`STRU`__ by a title __ELECTRON_NUMBER__ and a following line of the electron number value.

As for path integral, if there's a need to set different coordinates of beads at the beginning, set the centroid coordinates in __`STRU`__, and put different coordinates of beads in __`BEADS`__. For this situation, don't forget to add parameter __initial_beads__ into __`INPUT`__. The unit of coordinates in __`BEADS`__ is Bohr.

## __(work path)__
For the _DP_ situation, the folder isn't necessary. However, initial structure can also be put into a folder, and add the name of the folder with the parameter __work_path__ in __`INPUT`__.

For the _VASP_ situation, the folder is expected to contain input files of _VASP_, __`INCAR`, `KPOINTS`, `POTCAR`__ and any other necessary files, and the initial structure.

For the _ABACUS_ situation, the folder is expected to contain input files of _ABACUS_, __`INPUT`, `KPT`__ and any other necessary files, electron density file, __`SPIN1_CHG`__, and the initial structure.

## __(DP potential)__
For the _DP_ situation, don't forget to prepare a DP potential and refer to the potential with the parameter __DP_model__ in __`INPUT`__. However, for the _grand canonical_ ensemble, several DP potentials with different electron numbers should be prepared, put in the same potential, and named with the same prefix (refered by the parameter __DP_model_prefix__ in __`INPUT`__) following by the corresponding electron numbers.

# Output Files
### __`ALL_INPUT`__
The value of all parameters used in the simulation.

### __`energy.dat`__
Different physical quantities for each step.

### __`ALL_STRU`__
Dumped structures along the simulation process.

### __(`model_devi.out`)__
Model deviation results of virials and forces between different _DP_ potentials. Only output in the model deviation mode.

# Examples
### __MD simulation of C<sub>2</sub>H<sub>3</sub>__
MD simulation of C<sub>2</sub>H<sub>3</sub> under 79 K with a time interval 0.1 fs for 10000 steps.

_input files:_ [`INPUT`](./examples/MD/C2H3/INPUT), [`STRU`](./examples/MD/C2H3/STRU), `C2H3_compress.pb`

---
### __HMC simulation of C<sub>2</sub>H<sub>3</sub>__
HMC simulation of C<sub>2</sub>H<sub>3</sub> under 79 K for 200,000 steps.

_input files:_ [`INPUT`](./examples/HMC/C2H3/INPUT), `STRU`, `C2H3_compress.pb`

---
### __HMC simulation of C<sub>2</sub>H<sub>3</sub> in model deviation mode__
HMC simulation of C<sub>2</sub>H<sub>3</sub> under 79 K for 200,000 steps with model deviation on four different _DP_ potentials.

_input files:_ [`INPUT`](./examples/HMC/C2H3_model_devi/INPUT), `STRU`, `C2H3_compress_{0..3}.pb`

---
### __CHMC simulation of C<sub>2</sub>H<sub>3</sub>__
Constrained HMC simulation of C<sub>2</sub>H<sub>3</sub> with the "__difference__" reaction coordinate value 0.60 Bohr under 79 K for 200,000 steps.

_input files:_ [`INPUT`](./examples/CHMC/C2H3/INPUT), `STRU`, `C2H3_compress.pb`

---
### __CMC simulation of C<sub>2</sub>H<sub>3</sub>__
Constrained MC simulation of C<sub>2</sub>H<sub>3</sub> with the "__difference__" reaction coordinate value 0.60 Bohr under 79 K for 200,000 steps.

_input files:_ [`INPUT`](./examples/CMC/C2H3/INPUT), `STRU`, `C2H3_compress.pb`

---
### __CHMC simulation of monatomic molecule fluid system with LJ potential__
Constrained HMC simulation of monatomic molecule fluid system with the "__distance__" reaction coordinate value 0.95 Bohr under 300 K for 100,000 steps.

_input files:_ [`INPUT`](./examples/CHMC/LJ/INPUT), [`STRU`](./examples/CHMC/LJ/STRU)

> __Attention__: on account of the MC part realization for "__distance__" reaction coordinate simulation, the parameter __radial__ is required.

---
### __CHMC simulation of C<sub>3</sub>H<sub>4</sub>O<sub>2</sub>__
Constrained HMC simulation of C<sub>3</sub>H<sub>4</sub>O<sub>2</sub> via _ABACUS_ with the "__difference__" reaction coordinate value 0.60 Bohr under 300 K for 2,000 steps.

_input files:_ [`INPUT`](./examples/CHMC/C3H4O2_ABACUS/INPUT), `C3H4O2`/([`INPUT`](./examples/CHMC/C3H4O2_ABACUS/C3H4O2/INPUT), `KPT`, [`STRU`](./examples/CHMC/C3H4O2_ABACUS/C3H4O2/STRU), `OUT.ABACUS/SPIN1_CHG`)

> __Attention__: In the `INPUT` of _ABACUS_ `start_charge file` is necessary to use electron density of the previous single point calculation.

---
### __CHMC simulation of C<sub>3</sub>H<sub>4</sub>O<sub>2</sub>__
Constrained HMC simulation of C<sub>3</sub>H<sub>4</sub>O<sub>2</sub> via _VASP_ with the "__difference__" reaction coordinate value 0.69 Bohr under 300 K for 2,000 steps.

_input files:_ [`INPUT`](./examples/CHMC/C3H4O2_VASP/INPUT), `C3H4O2`/([`INCAR`](./examples/CHMC/C3H4O2_VASP/C3H4O2/INCAR), `KPOINTS`, `POTCAR`, `STRU`)

---
### __HMC simulation of Pt<sub>100</sub>-H<sub>25</sub>-(H<sub>2</sub>O)<sub>36</sub>__
HMC simulation of Pt<sub>100</sub>-H<sub>25</sub>-(H<sub>2</sub>O)<sub>36</sub> at 300 K for 100,000 steps under a _grand canonical_ ensemble condition. A DP potential involving the electron number is used. The electrochemical potential is set as -2.5 eV. A wall is set at the direction _x_ with coordinate 33 Bohr, and the masses of water are halved. 

_input files:_ [`INPUT`](./examples/GCHMC/PtHH2O/INPUT), [`STRU`](./examples/GCHMC/PtHH2O/STRU), `PtHH2O.pb`

> __Attention__: The electron number of the structure is put in the beginning of __`STRU`__ by a title __ELECTRON_NUMBER__ and a following line of the electron number value 0.54.

---
### __CHMC simulation of Pt<sub>100</sub>-H<sub>25</sub>-(H<sub>2</sub>O)<sub>36</sub>__
Constrained HMC simulation of Pt<sub>100</sub>-H<sub>25</sub>-(H<sub>2</sub>O)<sub>36</sub> with the "__difference__" reaction coordinate value -2.10 Bohr at 300 K for 100,000 steps under a _grand canonical_ ensemble condition. Five different DP potentials for the electron number ranging from 0.50 to 0.66 with the interval 0.04 is used. The electrochemical potential is set as -2.5 eV. A wall is set at the direction _x_ with coordinate 33 Bohr, and the masses of water are halved. 

_input files:_ [`INPUT`](./examples/GCCHMC/PtHH2O/INPUT), `STRU`, `compressed_0.50.pb`, `compressed_0.54.pb`, `compressed_0.58.pb`, `compressed_0.62.pb`, `compressed_0.66.pb`

---
### __CPIHMC simulation of C<sub>2</sub>H<sub>3</sub>__
Constrained PIHMC simulation of C<sub>2</sub>H<sub>3</sub> with the "__difference__" reaction coordinate value 1.20 Bohr under 79 K for 200,000 steps. The transfering proton is splited into 30 beads, and 2 beads changes in each PI step.

_input files:_ [`INPUT`](./examples/CHMC/C3H4O2_PI/INPUT), [`STRU`](./examples/CHMC/C3H4O2_PI/STRU), `C2H3_compress.pb`

If there's a need to set different coordinates of beads at the beginning, put different coordinates of beads in __`BEADS`__ and add parameter __initial_beads__ into __`INPUT`__.

_input files:_ [`INPUT`](./examples/CHMC/C3H4O2_PI/INPUT), [`STRU`](./examples/CHMC/C3H4O2_PI/STRU), [`BEADS`](./examples/CHMC/C3H4O2_PI/BEADS), `C2H3_compress.pb`

---
# Parameters
### __General Settings__
- __DP_model__
default: graph.pb (optional)

> the DP potential for the system. It's necessary when the DP potential is used for the potential energy and forces under a _canonical_ ensemble condition.

- __work_path__
default: . (optional)

> the work path of _ab initio_ calculation or the folder containing initial structure. Usually, it can be omitted when using the DP potential.

- __atom_file__
default: STRU (optional)

> the structure file. It's necessary only when you want to rename the structure file.

- __out_file__
default: energy.dat (optional)
                    
> the file to output different physical quantities.

- __calculation__
default: CHMC (optional)

> the simulation type.
> 
> choices: MD, HMC, CHMC, CMC

- __potential__
default: DP (optional)

> the style used to obtain the potential energy and forces.
> 
> choices: DP, VASP, ABACUS

- __OUTPUT__
default: ke pe etotal rc mfl mfr (optional)

> the physical quantities to be output in the file set by __out_file__. The physical quantities can be chosen contain kinetic energy, potential energy, total energy, internal energy, reaction coordinate, electron number, derivative of potential of mean force on the left atom, derivative of potential of mean force on the right atom, derivative of potential of mean force, accumulative average of potential energy, variance of potential energy.
>
> choices: ke pe etotal eint rc ne mfl mfr mf pe_ave pe_var

- __digits__
default: ke 6 pe 6 etotal 6 rc 2 mfl 6 mfr 6 (optional)

> the decimal digits of physical quantity(ies) to be output in the file set by __out_file__. The format of this parameter is string-integer pair(s) constituted by a physical quantity and a number defining the decimal digit of the physical quantity, such as _ke 2_.

- __steps__
(required)

> the number of simulation steps. A step means a sampling process between two adjacent acceptance judges.

- __dump__
default: 100 (optional)

> the number of steps in the interval between two outputs of the structures.

### __System Settings__
- __T__
default: 300

> unit: K
>                  
> the temperature of the simulation.

- __ntype__
(required)

> the number of element(s).

- __element_type__
(required)

> different element type(s) at the order in DP potential. The parameters are always required but only used for the DP potential when the parameter __element_index__ is not set.

- __element_index__
default: / (optional)

> different element index(s) used in DP potential mapping to element type(s) one-to-one. The parameter only displays when it has been set up.

- __RC_type__
default: Difference (optional)

> the reaction coordinate type.
>
> choices: Difference, Distance

- __reaction_coordinate__
default: 0

> unit: Bohr
> 
> the reaction coordinate value used in the constrained simulations. If not set manually, it will be set as the reaction coordinate value of the initial structure.

- __constraint_atom__
default: 0 0 0

> the atom indexes of the reaction coordinate. It's necessary for constrained simulations.
> 
> For the "__difference__" reaction coordinate, three atom indexes are required in the order $a_1$ $a_2$ $a_3$. Then the reaction coordinate value is $|\vec{r}_{a_2}-\vec{r}_{a_1}|-|\vec{r}_{a_2}-\vec{r}_{a_3}|$.
For the "__distance__" reaction coordinate, two atom indexes are required in the order $a_1$ $a_2$. Then the reaction coordinate value is $|\vec{r}_{a_2}-\vec{r}_{a_1}|$.

- __Wall__
default: None (optional)

> the hard wall(s) set to limit atom coordinates in the cell. The wall can be set on 6 directions. However, the function is not available in the MD simulations. 
>
> choices: xlo, xhi, ylo, yhi, zlo, zhi

- __W_pos__
default: / (optional)

> unit: Bohr
> 
> the position(s) of the set up wall(s). The number of the given number(s) should be equal to the number of the set up wall(s), and each number is mapped to each wall one-to-one. The parameter only works when the parameter __Wall__ has been set up.

- __N_MAX__
default: 1293 (optional)

> maximum number of atoms for the simulations, a temporarily useless parameter.

### __Grand Canonical Ensemble Settings__
- __DP_model_prefix__
default: NULL

> the prefix of DP potentials used for the _grand canonical_ ensemble. The simulation will be under a _grand canonical_ ensemble only when this parameter is set up.

- __GC_type__
default: Continuous

> the type of the DP potential's electron number.
>
> choices: Continuous, Discrete

- __Ne_range__
default: 0.5 0.66

> unit: 1
>
> the upper and lower limits of the electron number.

- __delta_Ne__
default: 0.04

> unit: 1
> 
> the change interval of the electron number.

- __Ne_proportion__
default: 0.1

> the proportion of the electron number change part in all moves.

- __mu__
default: -0.5

> unit: eV
> 
> the electrochemical potential for the _grand canonical_ ensemble.

### __Model Deviation Settings__
- __model_devi_models__
default: NULL

> the DP potentials used for the model deviations. More than one DP potentials should be provided, and the DP potential used for the simulation isn't contained inherently.

- __model_devi_file__
default: model_devi.out (optional)

> the file to output the model deviation results.

- __interval__
default: 100
> the interval between two model deviation results.

### __Molecular Dynamics Settings__
- __Delta_t__
default: 2.4

> unit: fs
>
> the time interval of Molecular Dynamics simulations.

### __Hybrid Monte Carlo Settings__
- __M_scaling__
default: 1 (optional)

> the mass scaling coefficient used in the HMC scheme. Decreasing the masses of all elements owns the similar effect of increasing the velocity sampling temperature.

- __Delta_t__
default: 2.4

> unit: fs
> 
> the time interval of Hybrid Monte Carlo simulations.

- __HMC_proportion__
default: 0.5

> the proportion of the HMC part in centroid move. The centroid move contains HMC part and MC part.

- __N_step__
default: 3 (optional)

> the number of V-V steps between two judges for whether accept or not.

- __HMC_lstep__
default: 0.3

> unit: Bohr
>
> the length threshold of MC movement in HMC method.

- __HMC_theta_step__
default: 0.1309

> unit: rad
>
> the theta threshold of MC movement in HMC method.

- __HMC_phi_step__
default: 1.0472

> unit: rad
> 
> the phi threshold of MC movement in HMC method.

### __Monte Carlo Settings__
- __radial__
default: 0.4

> unit: Bohr
>
> the radial of ball to adjust angle.

- __MC_theta_step__
default: 0.1309

> unit: rad
> 
> the theta threshold of central atom moving on the hyperboloid.

- __MC_phi_step__
default: 1.0472

> unit: rad
> 
> the phi threshold of central atom moving on the hyperboloid.

### __Path Integral Settings__
- __initial_beads__
default: (optional)

> the initial beads coordinates file. If the parameter isn't given, all beads are the same as the centroid at the beginning.

- __beads_file__
default: BEADS (optional)

> the prefix of beads coordinates file for the output structure.

- __P__
default: 1 (optional)

> the number of beads.

- __J__
default: 0 (optional)

> the number of beads changed in each step.

- __N_bead__
default: 1 (optional)

> the number of atom(s) to be considered with path integral.

- __bead_index__
default: 0 (optional)

> the index(s) of different bead atom(s) in coordinates.
