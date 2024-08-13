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

#ifndef DEM_PARTICLE_CLOSEST_POINT_H
#define DEM_PARTICLE_CLOSEST_POINT_H

inline Vector3d DEM_PARTICLE::GetClosestPoint2Sphere(Vector3d x)
{
	Vector3d cp = R*(x-X).normalized()+X;
	return cp;
}

inline Vector3d DEM_PARTICLE::GetClosestPoint2Cylinder(Vector3d x)
{
	double h0 = P0[P0.size()-2](2);
	double r0 = P0[0](0);
    Vector3d xr = Move2BodyFrame(x);
    Vector3d cpr = ClosestPoint::CylinderClosestPoint(xr, h0, r0);
    Vector3d cp = Move2GlobalFrame(cpr);
    return cp;
}

inline Vector3d DEM_PARTICLE::GetClosestPoint2Cuboid(Vector3d x)
{
    Vector3d xr = Move2BodyFrame(x);
    Vector3d cpr = ClosestPoint::CuboidClosestPoint(xr, P0[6]);
    Vector3d cp = Move2GlobalFrame(cpr);
    return cp;
}

// inline Vector3d DEM_PARTICLE::GetClosestPoint2Polygon2D(Vector3d x)
// {
// 	Vector3d xr = Qfi._transformVector(x - X);
//     return ClosestPoint::Polygon2DClosestPointToPoint(xr, P0, Edges);
// }

inline Vector3d DEM_PARTICLE::GetClosestPoint(Vector3d x)
{
	Vector3d cp;
	switch (ShapeType)
	{
	case 1 :
		cp = GetClosestPoint2Sphere(x);
		break;
	case 2 :	
		cp = GetClosestPoint2Cuboid(x);
		break;
	case 3 :	
		cp = GetClosestPoint2Cylinder(x);
		break;
	// case 7 :	
	// 	cp = GetClosestPoint2Polygon2D(x);
	// 	break;
	default :
		cout << "GetClosestPoint dont support this ShapeType yet" << endl;
		abort();
	}
	return cp;
}

#endif