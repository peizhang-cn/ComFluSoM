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

#include <DEM.h>

int main(int argc, char const *argv[])
{
	int nx = 80;
	int ny = 80;
	int nz = 100;
	DEM* a = new DEM(nx, ny, nz);

	a->Periodic[0] = true;
	a->Periodic[1] = false;
	a->Periodic[2] = false;
	
	a->Init();
	a->Nproc = 20;
	cout << "after init" << endl;
	cout << a->Flag[33][20][34] << endl;

	Vector3d x0 (0.,0.,0.);
	Vector3d x1 (nx,ny,nz);

	// Vector3d x0 (25.,25.,10.8);
	// Vector3d x1 (25.,25.,58.);

	Vector3d v0 (0.,0.00,-0.001);

	// Vector3d g (5.e-7,0.,-1.e-6);
	Vector3d g (0.,0.,-2.e-6);
	// Vector3d v1 (-0.00001,0.,0.);

	// Vector3d w0 (0.,0.,0.*M_PI/10000.);
	// Vector3d w1 (0.,0.,0.*M_PI/10000.);

	a->AddNSpheres(0, 700, x0, x1, 4, 0., 1.);

	// a->AddSphere(0, 10., x0, 1.0);
	// int ind = a->Lp.size()-1;
	// a->Lp[ind]->V = v0;
	// a->AddSphere(1, 10., x1, 1.0);
	// ind = a->Lp.size()-1;
	// a->Lp[ind]->Fix();
	// x0 << 31.6570100974, 12.5253798522, 26.2312868606;
	// a->AddSphere(0, 10., x0, 1.0);

	// x0 << 21.7211296385, 11.4343053682, 67.5540429394;
	// a->AddSphere(-1, 10., x0, 1.0);

	// x0 << 26.9376445804, 35.4862940276, 60.9602900473;
	// a->AddSphere(-2, 10., x0, 1.0);

	// x0 << 13.3110475812, 27.7388342716, 38.2121350892;
	// a->AddSphere(-3, 10., x0, 1.0);

	// x0 << 19.5437125045, 31.5124479758, 16.7693460927;
	// a->AddSphere(-4, 10., x0, 1.0);



	// x0 << 36.0872263764, 26.63474817, 41.9300870173;
	// a->AddSphere(-5, 10., x0, 1.0);

	// x0 << 11.9852524404, 11.7579841953, 16.8039474952;
	// a->AddSphere(-6, 10., x0, 1.0);


	a->SetG(g);
	// int cout = 1;
	// for (int i=1; i<(nx/2-1); ++i)
	// for (int j=1; j<(ny/2-1); ++j)
	// {
	// 	Vector3d x (2*i+1,2*j+1,1.5);
	// 	if (cout <5000000)
	// 	{
	// 		a->AddSphere(0, 0.999999, x, 1.);
	// 		ind = a->Lp.size()-1;
	// 		a->Lp[ind]->V = v0;
	// 		// a->Lp[cout]->Fix();
	// 		cout++;			
	// 	}
	// 	// a->AddSphere(0, 0.999999, x, 1.);
	// 	// ind = a->Lp.size()-1;
	// 	// a->Lp[ind]->V = v0;
	// 	// // a->Lp[cout]->Fix();
	// 	// cout++;
	// 	// x << 2*i+1,2*j+1,nz-1;
	// 	// a->AddSphere(1, 0.99, x, 1.0);
	// }

	// cout << "np= " << a->Lp.size() << endl;


	// a->Lp[0]->V = v0;
	// a->Lp[1]->V = v1;

	// // a->Lp[0]->Gn = 0.;
	// // a->Lp[1]->Gn = 0.;

	// // cout << a->Lrvt[0].i << endl;
	// // cout << a->Lrvt[0].j << endl;
	// // cout << a->Lrvt[0].Rvt << endl;
	// // abort();

	// a->Lp[0]->W = w0;
	// a->Lp[1]->W = w1;

	// Vector3d w (0., 0., 0.5*M_PI);
	// Vector3d r (1., 0., 0.);

	// cout << w.cross(r).transpose() <<endl;

	a->Solve(1.0e9, 20000, 0.02);

	return 0;
}