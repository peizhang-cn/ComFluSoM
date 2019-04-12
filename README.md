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
- Install dependencies (eigen3 and hdf5).
```
sudo apt install libeigen3-dev
sudo apt install libhdf5-dev
```
- Set environment variable.
```
sudo gedit ~/.bashrc
// Add following line to .bashrc
export ComFluSoM=~/ComFluSoM
```
# How to compile and run
- use make to compile and ./ to run.
