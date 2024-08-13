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

#ifndef DEM_PARTICLE_OBB_H
#define DEM_PARTICLE_OBB_H

inline void DEM_PARTICLE::InitOBB(size_t d, double e)
{
	Vector3d maxX = P0[0];
	Vector3d minX = P0[0];
	
	if (ShapeType==1)
	{
		OBB[0].setZero();
		OBB[1] << R+e, R+e, R+e;
	}
	else
	{
		for (size_t p=0; p<P0.size(); ++p)
		{
			if (P0[p](0)<minX(0))	minX(0) = P0[p](0);
			if (P0[p](1)<minX(1))	minX(1) = P0[p](1);
			if (P0[p](2)<minX(2))	minX(2) = P0[p](2);

			if (P0[p](0)>maxX(0))	maxX(0) = P0[p](0);
			if (P0[p](1)>maxX(1))	maxX(1) = P0[p](1);
			if (P0[p](2)>maxX(2))	maxX(2) = P0[p](2);
		}

		OBB[0] = 0.5*(minX+maxX);	// OBB centre
		OBB[1] << e, e, e;
		OBB[1] += 0.5*(maxX-minX);	// OBB half length
	}
	if (d==2)	OBB[1](2) = 0.;
}

#endif