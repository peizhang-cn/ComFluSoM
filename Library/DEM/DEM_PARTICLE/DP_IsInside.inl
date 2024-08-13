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

#ifndef DEM_PARTICLE_ISINSIDE_H
#define DEM_PARTICLE_ISINSIDE_H

inline bool DEM_PARTICLE::IsInsideSphere(Vector3d x)
{
	bool inside = true;
    if ((x-X).norm()>R)		inside = false;
    return inside;
}

inline bool DEM_PARTICLE::IsInsideCylinder(Vector3d x)
{
	bool inside = false;
	double h0 = P0[P0.size()-2](2);
	double r2 = P0[0](0)*P0[0](0);
	Vector3d xr = Qfi._transformVector(x - X);
	double hr = xr(2);
	double r2r = xr(0)*xr(0) + xr(1)*xr(1);
	if (abs(hr)<=h0 && r2r<=r2)	inside = true;
	return inside;
}

inline bool DEM_PARTICLE::IsInsideCuboid(Vector3d x)
{
    bool inside = true;
    Vector3d xr = Qfi._transformVector(x - X);

    double disx = abs(xr(0))-P0[6](0);
    if (disx>0.)    inside = false;
    else
    {
        double disy = abs(xr(1))-P0[6](1);
        if (disy>0.)    inside = false;
        else
        {
            double disz = abs(xr(2))-P0[6](2);
            if (disz>0.)    inside = false;
        }
    }
    return inside;
}

// inline bool DEM_PARTICLE::IsInsidePolygon2D(Vector3d x)
// {
// 	Vector3d xr = Qfi._transformVector(x - X);
//     return PointIsInsidePolygon2D(xr, P0, Edges);
// }

inline bool DEM_PARTICLE::IsInside(Vector3d x)
{
	bool inside = false;
	switch (ShapeType)
	{
	case 1 :
		inside = IsInsideSphere(x);
		break;
	case 2 :	
		inside = IsInsideCuboid(x);
		break;
	case 3 :	
		inside = IsInsideCylinder(x);
		break;
	// case 7 :	
	// 	inside = IsInsidePolygon2D(x);
	// 	break;
	default :
		cout << "CalDistance dont support this ShapeType yet" << endl;
		abort();
	}
	return inside;
}

#endif