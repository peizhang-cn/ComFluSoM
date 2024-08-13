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
 
#ifndef DEM_PARTICLE_DISTANCE_H
#define DEM_PARTICLE_DISTANCE_H

inline double DEM_PARTICLE::GetSignedDistance2Sphere(Vector3d x)
{
    double dis = (x-X).norm()-R;
    return dis;
}

inline double DEM_PARTICLE::GetSignedDistance2Cuboid(Vector3d x)
{
    Vector3d xr = Move2BodyFrame(x);
    double sdis = SignedDistance::CuboidSignedDistance(xr, P0[6]);
    return sdis;
}

inline double DEM_PARTICLE::GetSignedDistance2Cylinder(Vector3d x)
{
	double h0 = P0[P0.size()-2](2);
	double r0 = P0[0](0);
    Vector3d xr = Move2BodyFrame(x);
    double sdis = SignedDistance::CylinderSignedDistance(xr, h0, r0);
    return sdis;
}

// inline double DEM_PARTICLE::GetDistance2Polygon2D(Vector3d x)
// {
//     Vector3d xr = Qfi._transformVector(x - X);
//     double dis = Polygon2DDistance(xr, P0, Edges);
//     bool xIsInside = PointIsInsidePolygon2D(xr, P0, Edges);
//     if (xIsInside) dis *= -1.;
//     return dis;
// }

// efficiency is not great for this function
inline double DEM_PARTICLE::GetSignedDistance(Vector3d x)
{
	double dis = -1.;
	switch (ShapeType)
	{
	case 1 :
		dis = GetSignedDistance2Sphere(x);
		break;
	case 2 :	
		dis = GetSignedDistance2Cuboid(x);
		break;
	case 3 :	
		dis = GetSignedDistance2Cylinder(x);
		break;
	// case 7 :	
	// 	dis = GetDistance2Polygon2D(x);
	// 	break;
	default :
		cout << "GetDistance dont support this ShapeType yet" << endl;
		abort();
	}
	return dis;
}

inline double DEM_PARTICLE::GetDistance(Vector3d x)
{
	double sdis = GetSignedDistance(x);
	return abs(sdis);
}

#endif