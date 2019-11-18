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
	int nx = 200;
	int ny = 100;
	int nz = 0;
	// Create MPM domain
	MPM* a = new MPM(/*shape function type*/3, nx, ny, nz, gridSize);
	// Initialization
	a->Init();
	// Physcial parameters of particles
	double rhofPhysical 	= 1000.;			// Physical density, unit [kg/m^3]
	double rhosPhysical 	= 2000.;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force
	double miuPhysical 		= 1.0e-3;			// viscosity, unit [kg/(m*s)]
	double CsPhysical  		= 50.;				// Speed of sound, unit [m/s]
	double YoungPhysical 	= 1.0e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PossionPhysical  = 0.3;				// Possion ratio
	// Space time and mass step
	double dx = 0.02;							// unit [m]
	double dt = 2.0e-5;							// unit [s]
	double dm = 1.0e-1;							// unit [kg]
	// How many particles in a cell
	double Ratio = 1./4.;
	// Dementionless parameters of particles
	Vector3d G 		= GPhysical*dt*dt/dx;
	double Rhof 	= rhofPhysical*pow(dx,3)/dm;
	double Mpf 		= rhofPhysical*pow(dx,3)*pow(Ratio,2)/dm;
	double Mps 		= rhosPhysical*pow(dx,3)*pow(Ratio,2)/dm;
	double Miu 		= miuPhysical*dx*dt/dm;
	double Cs 		= CsPhysical*dt/dx;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Possion 	= PossionPhysical;
	// Start point of the box for generating particles
	Vector3d x0 (10, 10, 0);
	// demention of the box
	Vector3d l0 (30, 60, 0);
	a->Nproc = 30;
	a->Dc = 0.;
	// Define speed of sound
	a->Cs = Cs;
	// Generate a box of particles
	a->AddBoxParticles(-1, x0, l0, Ratio, Mpf);
	// Define gravity
	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->SetNewtonian(Miu);
		a->Lp[p]->B  = G;
		a->Lp[p]->P  = -Rhof*G(1)*(a->Lp[p]->X(1)-x0(1));
		a->Lp[p]->Stress(0,0) = a->Lp[p]->Stress(1,1) = a->Lp[p]->P;
	}
	size_t np0 = a->Lp.size();
	// Start point of the box for generating particles
	Vector3d x1 (40, 11, 0);
	// demention of the box
	Vector3d l1 (2, 25, 0);	
	a->AddBoxParticles(0, x1, l1, Ratio, Mps);
	for (size_t p=np0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->SetElastic(Young, Possion);
		a->Lp[p]->B  = G;
		// if (a->Lp[p]->X(1)>x1(1)+20.)	a->Lp[p]->FixV = true;
	}
	// Define left wall boundary
	for (int i=9; i<=10; ++i)
	for (int j=0; j<=ny; ++j)
	{
		Vector3d norm (-1., 0., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Define boundary
	for (int i=40; i<=42; ++i)
	for (int j=30; j<=35; ++j)
	{
		Vector3d norm (1., 0., 0.);
		a->SetNonSlippingBC(i,j,0);
	}
	// Define middle half wall boundary
	for (int i=40; i<=41; ++i)
	for (int j=30; j<=75; ++j)
	{
		Vector3d norm (1., 0., 0.);
		Vector3d norm1 (-1., 0., 0.);
		a->SetSlippingBC(i,j,0, norm);
		a->SetSlippingBC(i,j,0, norm1);
	}
	// Define right wall boundary
	for (int i=149; i<=150; ++i)
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
	for (int j=99; j<=100; ++j)
	{
		Vector3d norm (0., 1., 0.);
		a->SetSlippingBC(i,j,0, norm);
	}
	// Solve
	a->SolveMUSL(/*total time step*/200000,/*save per x time step*/500);
	return 0;
}