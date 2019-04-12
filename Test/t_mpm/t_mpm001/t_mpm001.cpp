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
	Vector3d dd (1,1,1);

	int nx = 100;
	int ny = 100;
	int nz = 100;
	MPM* a = new MPM(3, nx, ny, nz, dd);

	a->Nproc = 1;

	a->Init();

	Vector3d x0 (40, 40, 18);
	Vector3d l0 (4, 4, 22);

	double rhosPhysical 	= 1000.;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., -9.8, 0.);			// Body force

	double YoungPhysical 	= 5.0e6;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.3;				// Poisson ratio

	double dx = 0.1;							// unit [m]
	double dt = 0.5e-3;							// unit [s]
	double dm = 1.0e-4;							// unit [kg]
	double Ratio = 1./3.;

	Vector3d G 		= GPhysical*dt*dt/dx;
	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,3)/dm;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Poisson 	= PoissonPhysical;

	a->AddBoxParticles(x0, l0, Ratio, Mp, Young, Poisson);

	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->B  = G;
	}

	for (int i=0; i<nx; ++i)
	for (int j=0; j<ny; ++j)
	for (int k=0; k<=20; ++k)
	{
		Vector3i cell {i,j,k};
		a->LFn.push_back(cell);
	}

	a->SolveMUSL(100000,100);
	// a->SolveUSF(10000,100);

	return 0;
}