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

// 2D simlation of elastic beam under gravity

#include <MPM.h>

int main(int argc, char const *argv[])
{
	// Size of background grid
	double dx = 0.02;					// space step [m]
	double dt = 1.e-4;					// time step [s]
	Vector3d origin (0.,0.,0.);			// origin of the background grid [m]
	Vector3d lx (4.,4.,0.);				// size of the background grid [m]

	// Physcial parameters of particles
	double rho 		= 3000.;			// Physical density, unit [kg/m^3]
	double young 	= 1.5e6;			// Young’s modulus [Pa]
	double poisson	= 0.3;				// Poisson ratio
	Vector3d g (0., -9.8, 0.);			// Body force

	// Create MPM domain
	const int SType = 1;
	const int D = 2;
	MPM<SType, D>* mpm = new MPM<SType, D>(origin, lx, dx, dt);

	mpm->Dc = 0.2;

	// Initialization
	bool useFbar = false;
	mpm->Init(useFbar);

	// choose how many processors to use
	mpm->SetParallel(1);

	// Start point of the box for generating particles
	Vector3d x0 (1.1, 1.1, 0.);
	// Ending point of the box
	Vector3d x1 (1.5, 1.2, 0.);
	mpm->AddBoxParticles(-1, x0, x1, 0.01, rho);

	// Set viscosity and gravity
	for (size_t p=0; p<mpm->Lp.size(); ++p)
	{
		mpm->Lp[p]->SetLinearElastic(young, poisson);
		mpm->Lp[p]->B  = g;
	}

	for (size_t n=0; n<mpm->Ln.size(); ++n)
	{
		// left wall
		if (mpm->Ln[n]->X(0)<=x0(0))
		{
			mpm->SetNonSlippingBC(n);
		}
	}

	mpm->SolveMUSL(/*total time step*/80000,/*save per x time step*/200);
	return 0;
}