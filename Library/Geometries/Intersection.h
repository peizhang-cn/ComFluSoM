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

#ifndef INTERSECTION_H
#define INTERSECTION_H

// this function finds line line interscetion u (s,t), where the interecion point is A+s*(B-A), C+t*(D-C)
template <typename T>
Vector2d LineLineIntersction2D(T A, T B, T C, T D)
{
	T AB = B-A;
	T DC = C-D;
	T AC = C-A;
	Matrix2d M;
	// M(0,0) = CD(0); M(1,0) = CD(1);
	// M(0,1) = -AB(0); M(1,1) = -AB(1);

	// Vector2d b (CA(0), CA(1));
	M(0,0) = AB(0); M(1,0) = AB(1);
	M(0,1) = DC(0); M(1,1) = DC(1);

	Vector2d b (AC(0), AC(1));
	// Vector2d u = M.colPivHouseholderQr().solve(b);
	Vector2d u = M.inverse()*b;
	return u;
}

// this function finds line triangle interscetion by using Barycentric Coordinates
Vector3d LineTriangleIntersction(Vector3d x0, Vector3d x1, Vector3d p0, Vector3d p1, Vector3d p2)
{
	Matrix3d A;
	A.col(0) = p0-p2;
	A.col(1) = p1-p2;
	A.col(2) = x0-x1;

	Vector3d u (-1.,-1.,-1.);
	Vector3d n = A.col(0).cross(A.col(1));
	if (n.dot(A.col(2))!=0.)
	{
		Vector3d b = x0-p2;
		// u = A.colPivHouseholderQr().solve(b);
		u = A.inverse()*b;
	}
	// else
	// {
	// 	cout << "wrong in LineTriangleIntersction" << endl;
	// 	// abort();
	// }
	return u;
}

#endif