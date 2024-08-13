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

#ifndef AREA_H
#define AREA_H

// calculate the area of a non-self-intersect polygon in 2d 
template <typename T1, typename T2>
double Polygon2DArea(vector<T1> P, vector<T2> E)
{
	double area = 0.;
	for (size_t e=0; e<E.size(); ++e)
	{
		T1 A = P[E[e](0)];
		T1 B = P[E[e](1)];
		area += -A(1)*B(0) + A(0)*B(1);
	}
	return 0.5*abs(area);
}

#endif