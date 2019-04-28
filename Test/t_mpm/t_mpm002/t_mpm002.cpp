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
	MPM* a = new MPM(0, 1, nx, ny, nz, dd);

	a->Nproc = 1;

	a->Init();

	Vector3d x0 (40, 40, 20);
	Vector3d l0 (20, 20, 20);

	Vector3d x1 (10, 10, 19);
	Vector3d l1 (80, 80, 1);

	// Vector3d ratio (1./3., 1./3., 1.);

	double rhosPhysical 	= 2039.435;			// Physical density, unit [kg/m^3]
	Vector3d GPhysical (0., 0., -9.8);			// Body force

	double YoungPhysical 	= 7.5e7;			// Young's modus, unit [kg/(m*s^2)] (or Pa)
	double PossionPhysical  = 0.3;				// Possion ratio

	double CPhysical		= 5.0e3;			// Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double PhiPhysical		= 30./180.*M_PI;	// Angle of internal friction
	double PsiPhysical		= 0./180.*M_PI;		// Angle of dilatation

	double dx = 2.0;							// unit [m]
	double dt = 1.0e-3;							// unit [s]
	double dm = 1.0e-6;							// unit [kg]
	double Ratio = 1./3.;

	Vector3d G 		= GPhysical*dt*dt/dx;
	double Mp 		= rhosPhysical*pow(dx,3)*pow(Ratio,3)/dm;
	double Young 	= YoungPhysical*dx*dt*dt/dm;
	double Possion 	= PossionPhysical;

	double Rhos 	= rhosPhysical/dm*dx*dx*dx;

	double C 		= CPhysical*dx*dt*dt/dm;
	double Phi 		= PhiPhysical;
	double Psi 		= PsiPhysical;

	double K0 		= Possion/(1.-Possion);

	cout << "Mp= " << Mp << endl;
	cout << "Young= " << Young << endl;
	cout << "Possion= " << Possion << endl;
	cout << "G= " << G.transpose() << endl;
	cout << "C= " << C << endl;
	// double Young 	= 1.0e2;
	// double Possion 	= 0.4;
	// double Mu = 0.5*Young/(1.+Possion);
	// double La = Young*Possion/(1+Possion)/(1-2*Possion);

	a->AddBoxParticles(x0, l0, Ratio, Mp, Young, Possion);

	Vector3d xc (50., 50., 0);

    Vector3d x_sph = x0+0.5*l0;

    // x_sph(0) += 2.;
    // x_sph(2) += 2.;

    double r0 = 5.;
    double r1 = 42.;

	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		double h = -a->Lp[p]->X(0) + 80.;

		if (a->Lp[p]->X(2)>h && (a->Lp[p]->X-x_sph).norm()>r0)	a->Lp[p]->Removed = true;
		// else								a->Lp[p]->Removed = false;
	}
	a->DeleteParticles();

/*	a->WriteFileH5(0);
	a->WriteFileH5(1);

	abort();*/

	// a->AddBoxParticles(x1, l1, Ratio, Mp, Young, Possion);

	for (size_t p=0; p<a->Lp.size(); ++p)
	{
		a->Lp[p]->C = C;
		a->Lp[p]->Phi = Phi;
		a->Lp[p]->Psi = Psi;

		// a->Lp[p]->S(2,2) = (a->Lp[p]->X(2)-(x0(2)+l0(2)))*G(2)*Rhos;
		// a->Lp[p]->S(1,1) = a->Lp[p]->S(2,2)*K0;
		// a->Lp[p]->S(0,0) = a->Lp[p]->S(1,1);

		a->Lp[p]->B  = G;
		a->Lp[p]->FixV = true;

		if ((a->Lp[p]->X-x_sph).norm()<r0)
		{
			a->Lp[p]->Tag = -2;
			// if ((a->Lp[p]->X-x_sph).norm()>r1)
			// {
			// 	a->Lp[p]->Tag = -3;
			// 	a->Lp[p]->C = 0.;
			// 	// cout << "c= " << a->Lp[p]->C << endl;
			// }

			// a->Lp[p]->Tag = -3;
			a->Lp[p]->C = 0.0001*C;
			a->Lp[p]->FixV = false;
		}
		else
		{
			a->Lp[p]->Tag = -1;
		}
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