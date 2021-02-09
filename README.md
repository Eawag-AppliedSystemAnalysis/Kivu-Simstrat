# Kivu-Simstrat: a one-dimensional physical lake model including a biochemical library, optimized for the use on Lake Kivu

Simstrat is a one-dimensional physical lake model for the simulation of stratification and mixing in deep stratified lakes. The model was originally developed by Goudsmit et al. (2002) and has been successfully applied to lakes with different physical properties. A k-ε model is used to model turbulent mixing including energy transfer of internal seiches. River or groundwater inflow can be added at specific depths or as density-dependent intrusions. The newest version of Simstrat (see below) can also simulate ice/snow covers.

AED2 is a biogeochemical library developed by the University of Western Australia (https://github.com/AquaticEcoDynamics/libaed2)

Kivu-Simstrat consists of the coupled version of Simstrat-AED2, with some modifications for Lake Kivu.

The Simstrat documentation is in the `doc` folder. The AED2 documentation is pending, but the configuration file (aed2_config.nml) is well documented.
 
## Run Simstrat (only windows)
Pre-built windows binaries are available in the `build` folder.

To run the scenarios included in the Lake Kivu modelling publication, go to the root (Kivu_Simstrat_v1.0) and run

~~~bash
.\run_kivu_steady.bat
~~~

for a 500 year steady-state simulation, starting with T and S from Ross et al., 2015 (PlosOne) and CH4 and CO2 from Bärenbold et al., 2020 (PlosOne).
The results of this run are stored in scenarios/Results_steady. The units of results are the same as in inflow/initial conditions.


Run

~~~bash
.\run_kivu_steady_ch4inflow.bat
~~~

for a 500 year steady-state simulation, starting with T and S from Ross et al., 2015 (PlosOne) and CH4 and CO2 from Bärenbold et al., 2020 (PlosOne), with 50% of CH4 from groundwater inflow.
The results of this run are stored in scenarios/Results_steady. The units of results are the same as in inflow/initial conditions.


Run

~~~bash
.\run_kivu_longterm.bat
~~~

for a 2000 year longterm simulation, starting with a homogeneous lake.
The results of this run are stored in scenarios/Results_longterm. The units of results are the same as in inflow/initial conditions.


Run

~~~bash
.\run_kivu_longterm_ch4_inflow.bat
~~~

for a 2000 year longterm simulation, starting with a homogeneous lake and with 50 % of CH4 coming from inflow.
The results of this run are stored in scenarios/Results_longterm. The units of results are the same as in inflow/initial conditions.


Run

~~~bash
.\run_kivu_tritium.bat
~~~

for a 65 year Tritium (3H) simulation, starting with 3H=0 and forced by RCP6 timeseries. The 3H data is published in a separate data package.
The results of this run are stored in scenarios/Results_tritium. The units of results are the same as in inflow/initial conditions.




## Build AED2 and Simstrat yourself (only windows)

Please install the required packages listed below (System requirements) and then from terminal first navigate into `lib\libaed2` folder and run

~~~bash
mingw32-make
~~~

Then, move to the `build` folder and run 

~~~bash
FoBiS.py build
~~~

**System requirements**

- [Python](https://www.python.org/) 2.7 or later
- [FoBiS.py](https://github.com/szaghi/FoBiS) 2.2.6 or later (available via GitHub or pip using `pip install FoBiS.py`)
- 2 compiler options:
	- [Intel Fortran (Intel Parallel Studio XE 2016)](https://software.intel.com/en-us/parallel-studio-xe/choose-download) (commercial)
	- Gfortran 6.3 or later (free)