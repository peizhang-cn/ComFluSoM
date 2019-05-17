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
	int nz = 0;
	MPM* a = new MPM(3, 2, nx, ny, nz, dd);

	a->Nproc = 1;

	a->Init();

	Vector3d x0 (5, 5, 0);
	Vector3d l0 (50, 25, 0);

	double rhosPhysical 	= 1000.;			// Physical density, unit [kg/m^3]
	double YoungPhysical 	= 3.0e5;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.;				// Poisson ratio

	double dx = 0.01;							// unit [m]
	double dt = 2.e-5;							// unit [s]
	double dm = 1.0e-0;							// unit [kg]
	double Ratio = 1./3.;

	// GPhysical(1) = 700./(Mp*dm*a->Lp.size());

	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,3)/dm;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Poisson 	= PoissonPhysical;

	// double Mp 		= 8000./288.;
	// double Young 	= 3.e5;
	// double Poisson 	= 0.;

	a->AddBoxParticles(x0, l0, Ratio, Mp, Young, Poisson);

	// Vector3d x1 (1, 1, 0);
	// Vector3d l1 (90, 5, 0);

	// a->AddBoxParticles(x1, l1, Ratio, Mp, Young, Poisson);

	Vector3d GPhysical (0., -9.8, 0.);
	Vector3d G 		= GPhysical*dt*dt/dx/dm;

	double MuPhysical 	= 1.0e-3; 				// Dynamic viscosity, unit [kg/(m*s)] (or Pa*s)
	double Mu 			= MuPhysical*dx*dt/dm;

	double CPhysical 	= 50.;
	double C 			= CPhysical*dt/dx;

	// Vector3d GPhysical (0., 0., 0.);
	// GPhysical(2) = -800.;
	// Vector3d G 		= GPhysical*dt*dt/dx;
	// G = GPhysical;

	cout << "Mp= " << Mp << endl;
	cout << "G= " << G.transpose() << endl;

	a->Dt = 1.;
	a->C = C;

	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->Mu  = Mu;
		a->Lp[p]->B  = G;

		// if (a->Lp[p]->X(1)<6.)		a->Lp[p]->FixV = true;
	}

	for (int i=0; i<=5 ; ++i)
	for (int j=0; j<=ny; ++j)
	{
		Vector3i cellL {i,j,0};
		a->LFn.push_back(cellL);
		Vector3i cellR {nx-i,j,0};
		a->LFn.push_back(cellR);
	}

	for (int i=0; i<=nx; ++i)
	for (int j=0; j<=5 ; ++j)
	{
		Vector3i cellB {i,j,0};
		a->LFn.push_back(cellB);
		Vector3i cellT {ny-i,j,0};
		a->LFn.push_back(cellT);
	}

	a->SolveMUSL(1000000,100);
	// a->SolveUSF(10000,100);

	return 0;
}