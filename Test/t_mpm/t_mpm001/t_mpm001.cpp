/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of one grid
	Vector3d gridSize (1,1,1);
	// Domain size
	int nx = 500;
	int ny = 500;
	int nz = 0;
	// Create MPM domain
	MPM* a = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// switch MLS for velocity interpolation
	a->MLSv = true;
	// Initialization
	a->Init();
	// Physcial parameters of particles
	double rhosPhysical 	= 2039.435;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double YoungPhysical 	= 7.5e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PossionPhysical  = 0.3;				// Possion ratio
	// Space time and mass step
	double dx = 0.5;							// unit [m]
	double dt = 1.0e-4;							// unit [s]
	double dm = 1.0e-1;							// unit [kg]
	// How many particles in a cell
	double Ratio = 1./4.;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,2)/dm;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Possion 	= PossionPhysical;
	// Start point of the box for generating particles
	Vector3d x0 (240, 240, 0);
	// demention of the box
	Vector3d l0 (20, 4, 0);
	a->Nproc = 1;
	a->Dc = 0.;
	// Generate a box of particles
	// a->AddBoxParticles(x0, l0, Ratio, Mp, Young, Possion);
	a->AddBoxParticles(-1, x0, l0, Ratio, Mp);
	// Define gravity
	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->SetElastic(Young, Possion);
		a->Lp[p]->B  = G;
	}
	// Define boundary
	for (int i=239; i<=240; ++i)
	for (int j=0; j<ny; ++j)
	{
		Vector3d norm (1., 0., 0.);
		a->SetNonSlippingBC(i,j,0);
	}
	// Solve
	a->SolveMUSL(/*total time step*/50000,/*save per x time step*/100);
	return 0;
}