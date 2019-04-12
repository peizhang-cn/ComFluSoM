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
	DEM* a = new DEM(100, 100, 100);

	a->Init();

	Vector3d x0 (39.,50.,50.);
	Vector3d x1 (61.,50.,50.);

	Vector3d v0 (0.001,0.,0.);
	Vector3d v1 (-0.001,0.,0.);

	Vector3d w0 (0.,0.,2.*M_PI/10000.);
	Vector3d w1 (0.,0.,-2.*M_PI/10000.);

	a->AddSphere(0, 10., x0, 1.0);
	a->AddSphere(1, 10., x1, 1.0);

	a->Lp[0]->V = v0;
	a->Lp[1]->V = v1;

	a->Lp[0]->W = w0;
	a->Lp[1]->W = w1;

	a->Solve(10000, 100, 1.);
	return 0;
}