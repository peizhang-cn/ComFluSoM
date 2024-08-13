/****************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics                *
 * Copyright (C) 2024 Pei Zhang                                             *
 * Email: peizhang.hhu@gmail.com                                            *
 *                                                                          *
 * This program is free software: you can redistribute it and/or modify     *
 * it under the terms of the GNU Affero General Public License as           *
 * published by the Free Software Foundation, either version 3 of the       *
 * License, or (at your option) any later version.                          *
 *                                                                          *
 * This program is distributed in the hope that it will be useful,          *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU Affero General Public License for more details.                      *
 *                                                                          *
 * You should have received a copy of the GNU Affero General Public License *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.   *
 * In cases where the constraints of the Open Source license prevent you 	*
 * from using ComFluSoM, please contact by peizhang.hhu@gmail.com for a 	*
 * commercial license. 														*
 ****************************************************************************/

#ifndef DEM_CONTACT_SPHERE_CUBOID_H
#define DEM_CONTACT_SPHERE_CUBOID_H

inline void DEM::Sphere2Cuboid(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, vector<CONTACT_INFO>& lci)
{
	DEM_PARTICLE* pi = Pi;
	DEM_PARTICLE* pj = Pj;

	Vector3d xi = Xi;
	Vector3d xj = Xj;
	// make sure that pi is sphere and pj is cuboid
	double swapSign = 1.;
	if (pi->ShapeType!=1)
	{
		swap(pi,pj);
		swap(xi,xj);
		swapSign = -1.;
	}

	Vector3d cpj = pj->GetClosestPoint2Cuboid(xi);

	// cout << "cpj: " << cpj.transpose() << endl;

	Vector3d n = xi-cpj;							// Normal direction (pj pinnts to pi)
	double delta = pi->R+pj->Rs-n.norm(); 			// Overlapping distance for sphere
	
	// cout << "Delta: " << delta << endl;
	// abort();

	if (delta>0.)
	{
		CONTACT_INFO ci;
		ci.Delta = delta;
		n.normalize();
		ci.Xc = cpj+(pj->Rs-0.5*delta)*n;
		ci.Nc = swapSign*n;
		lci.push_back(ci);
	}
}

#endif