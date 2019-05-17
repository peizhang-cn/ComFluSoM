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
	MPM* a = new MPM(3, 0, nx, ny, nz, dd);

	a->Nproc = 1;

	a->Init();

	Vector3d x0 (50, 50, 50);
	Vector3d l0 (24, 3, 3);

	double rhosPhysical 	= 1000.;			// Physical density, unit [kg/m^3]
	double YoungPhysical 	= 3.0e5;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PoissonPhysical  = 0.;				// Poisson ratio

	double dx = 0.5;							// unit [m]
	double dt = 3.e-3;							// unit [s]
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

	Vector3d FPhysical (0., 0., -200.);
	Vector3d F 		= FPhysical*dt*dt/dx/dm;

	// Vector3d GPhysical (0., 0., 0.);
	// GPhysical(2) = -800.;
	// Vector3d G 		= GPhysical*dt*dt/dx;
	// G = GPhysical;

	cout << "Mp= " << Mp << endl;
	cout << "F= " << F.transpose() << endl;

	double dy = pow(l0(0),3)*F(2)/(3.*Young*pow(l0(1),4)/12.);

	cout << "dy= " << dy << endl;
	a->Dt = 1.;
	a->Dc = 0.85;

	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		if (a->Lp[p]->X(0)>73.5 && a->Lp[p]->X(1)==51.5 && a->Lp[p]->X(2)==51.5)
		{
			a->Lp[p]->Fh  = F;
			cout << "set g" << endl;
			cout << a->Lp[p]->X.transpose() << endl;
			cout << "p= " << p << endl;
		}
	}

	for (int i=0; i<51; ++i)
	for (int j=0; j<ny; ++j)
	for (int k=0; k<nz; ++k)
	{
		Vector3i cell {i,j,k};
		a->LFn.push_back(cell);
	}

	a->SolveMUSL(1000000,1000);
	// a->SolveUSF(10000,100);

	return 0;
}