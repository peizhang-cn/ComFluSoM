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

// 2D simlation of Newtonian fluid

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of one grid
	Vector3d gridSize (1,1,1);
	// Domain size
	int nx = 200;
	int ny = 100;
	int nz = 0;
	// Create MPM domain
	MPM* a = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// Initialization
	a->Init();
	// Physcial parameters of particles
	double rhopPhysical 	= 3000.;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double CsPhysical  		= 50.;				// Speed of sound, unit [m/s]
	double Phi = 0.53;							// Initial volume fraction
	double DPhysical = 2.0e-4; 					// Diameter of particles [m]
	// Space time and mass step
	double dx = 0.02;							// unit [m]
	double dt = 2.0e-5;							// unit [s]
	double dm = 1.0e-1;							// unit [kg]
	// How many particles in a cell
	double Ratio = 1./3.;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Rhop 	= rhopPhysical*pow(dx,3)/dm;
	double Rhos 	= Phi*Rhop;
	double Mp 		= Rhos*pow(Ratio,2);
	double Cs 		= CsPhysical*dt/dx;
	double D 		= DPhysical/dx;
	// Start point of the box for generating particles
	Vector3d x0 (10, 10, 0);
	// demention of the box
	Vector3d l0 (10, 5, 0);
	a->Nproc = 1;
	a->Dc = 0.;
	// Define speed of sound
	a->Cs = Cs;
	// Generate a box of particles
	a->AddBoxParticles(-1, x0, l0, Ratio, Mp);
	// Define gravity
	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->SetGranular(Rhop, D);
		a->Lp[p]->B  = G;
		a->Lp[p]->P  = -Rhos*G(1)*(a->Lp[p]->X(1)-x0(1));
		a->Lp[p]->Stress(0,0) = a->Lp[p]->Stress(1,1) = a->Lp[p]->P;
	}
	// Define boundary
	for (int i=9; i<=10; ++i)
	for (int j=0; j<=ny; ++j)
	{
		Vector3d norm (-1., 0., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Define boundary
	for (int i=170; i<=171; ++i)
	for (int j=0; j<=ny; ++j)
	{
		Vector3d norm (1., 0., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Define boundary
	for (int i=0; i<=nx; ++i)
	for (int j=9; j<=10; ++j)
	{
		Vector3d norm (0., -1., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Define boundary
	for (int i=0; i<=nx; ++i)
	for (int j=50; j<=51; ++j)
	{
		Vector3d norm (0., 1., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Solve
	a->SolveMUSL(/*total time step*/3000,/*save per x time step*/100);
	return 0;
}