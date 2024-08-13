# ComFluSoM
ComFluSoM is a simulation kit for **Com**plex **Flu**id, **So**lid and **So**il **M**echanics problems. Typically, Fluid-Structure Interactions, Fluid-Particle interactions.
# Methods
- **Lattice Boltzmann Method (LBM)** for flow.
- **Material Point Method (MPM)** for solid and soil deformations.
- **Discrete Element Method (DEM)** for rigid particle motions. 
- **Random Walking Method (RWM)** for solute (heat) transports.
# Enveriment
- Ubuntu 18.04
# How to install
- Download the source code.
- Install dependencies (eigen3, hdf5 and cgal).
```
sudo apt install libeigen3-dev
sudo apt install libhdf5-dev
sudo apt-get install libcgal-dev
```
- Make sure you have cmake version later than 3.3
# How to compile and run
- compile
- go to Test folder and into build under demo cases
```
cmake ..
make
```
- use ./ to run, for instance:
```
./t_dem001
```
