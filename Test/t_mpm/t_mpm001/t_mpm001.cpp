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

// 2D simlation of Newtonian fluid, dam break in a closed box

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of background grid
	double dx = 0.02;					// space step [m]
	double dt = 1.e-5;					// time step [s]
	Vector3d origin (0.,0.,0.);			// origin of the background grid [m]
	Vector3d lx (4.,2.,0.);				// size of the background grid [m]

	// Physcial parameters of particles
	double rho 		= 1000.;			// Physical density, unit [kg/m^3]
	double miu 		= 1.01e-3;			// viscosity, unit [kg/(m*s)]
	double cs  		= 50.;				// Speed of sound, unit [m/s]
	Vector3d g (0., -9.8, 0.);			// Body force

	// Create MPM domain
	const int SType = 1;
	const int D = 2;
	MPM<SType, D>* mpm = new MPM<SType, D>(origin, lx, dx, dt);

	// Initialization
	bool useFbar = true;
	mpm->Init(useFbar);

	// Set speed of sound
	mpm->Cs = cs;
	// choose how many processors to use
	mpm->SetParallel(8);

	// Start point of the box for generating particles
	Vector3d x0 (0.2, 0.2, 0.);
	// Ending point of the box
	Vector3d x1 (1.4, 0.8, 0.);
	mpm->AddBoxParticles(-1, x0, x1, 0.01, rho);

	// Set viscosity and gravity
	for (size_t p=0; p<mpm->Lp.size(); ++p)
	{
		mpm->Lp[p]->SetNewtonian(miu);
		mpm->Lp[p]->B  = g;
	}

	for (size_t n=0; n<mpm->Ln.size(); ++n)
	{
		// bottom wall, make sure two layers of nodes are included
		if (mpm->Ln[n]->X(1)<=x0(1) && mpm->Ln[n]->X(1)>=0.18)
		{
			Vector3d norm (0., -1., 0.);
			mpm->SetSlippingBC(n, norm);
		}
		// to wall, make sure two layers of nodes are included
		if (mpm->Ln[n]->X(1)>=1. && mpm->Ln[n]->X(1)>=1.02)
		{
			Vector3d norm (0., 1., 0.);
			mpm->SetSlippingBC(n, norm);
		}
		// left wall
		if (mpm->Ln[n]->X(0)<=x0(0) && mpm->Ln[n]->X(0)>=0.18)
		{
			Vector3d norm (-1., 0., 0.);
			mpm->SetSlippingBC(n, norm);
		}
		// right wall
		if (mpm->Ln[n]->X(0)>=2.4 && mpm->Ln[n]->X(0)<=2.42)
		{
			Vector3d norm (1., 0., 0.);
			mpm->SetSlippingBC(n, norm);
		}
	}

	mpm->SolveMUSL(/*total time step*/800000,/*save per x time step*/2000);
	return 0;
}