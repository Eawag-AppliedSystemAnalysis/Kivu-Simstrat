## Kivu-Simstrat Version 1.0 (November 2020)
**Architecture**
- Object-oriented, modern Fortran 2003 architecture
- JSON formatted configuration files
- Consistent array indexing in fortran standard (arrays start with index 1)
- Docker container to build the model

**Model update**
- Based on Simstrat 2.2
- Coupled with AED2
- Modified for Lake Kivu:
	- parameterization of double diffusion using data from Sommer et al., 2013
	- including effect of CO2 and CH4 ond water density
	- Simulation of water tracers 2H, 3H, 18O and noble gases He, Ne, 36Ar, Kr