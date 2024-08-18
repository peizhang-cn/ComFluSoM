# ComFluSoM
ComFluSoM is a simulation kit for **Com**plex **Flu**id, **So**lid and **So**il **M**echanics problems. Typically, Fluid-Structure Interactions, Fluid-Particle interactions.
# Methods
- **Lattice Boltzmann Method (LBM)** for flow.
- **Material Point Method (MPM)** for solid and soil deformations.
- **Discrete Element Method (DEM)** for rigid particle motions. 
- **Random Walking Method (RWM)** for solute (heat) transports.
# Current status
- ComFluSoM is under a major modification, thus only part of DEM code (with limited capability) is aviable. the MPM part will be released soon. openMP is used to parallelize on CPU, a GPU version is also in long-term plan.
- **DEM**
- contact model:

| linear | hertz |
|----------|----------|

- contact type:

| sphere-sphere | sphere-cuboid |
|----------|----------|


# Future update
- **DEM**
- contact type:

| sphere-polyhedra | convex polyhedra-polyhedra | convex metaball-metaball |
|----------|----------|----------|

- boundary condition:

| complex terrain |
|----------|

- **MPM**
- shape function:

| GIMP | B-spline |
|----------|----------|

- stabilization:

| F-bar |
|----------|

- constitutive model:

| elastic | mohr-coulomb | drucker-prager | newtonian |
|----------|----------|----------|----------|

- boundary condition:

| no-slip | free-slip | friction | complex terrain |
|----------|----------|----------|----------|

- **LBM**
- collision model:

| SRT | MRT | Entropic |
|----------|----------|----------|
  
- boundary condition:

| no-slip | free-slip | velocity | pressure |
|----------|----------|----------|----------|

- **DEM-LBM**
- coupling scheme:

| IBB | IBM |
|----------|----------|

- **DEM-MPM**
- coupling scheme:

| momentum exchange |
|----------|
 
- boundary condition:

| no-slip | free-slip |
|----------|----------|

# Enveriment
- Tested on Ubuntu 22.04, latter version of Ubuntu should also work.
# How to install
- Download the source code.
- Install dependencies (eigen3, hdf5 and cgal).
```
sudo apt install libeigen3-dev
sudo apt install libhdf5-dev
sudo apt-get install libcgal-dev
```
- Make sure you have cmake version later than 3.3.
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
# Visualization
- Visit (https://visit-dav.github.io/visit-website/index.html)
- paraview should also work, but I don't use it personally 
# license
- Copyright (C) 2024 Pei Zhang
- Email: peizhang.hhu@gmail.com
- ComFluSoM is an open-source software under AGPL license, for commertial use, please contact me for a dual license.
