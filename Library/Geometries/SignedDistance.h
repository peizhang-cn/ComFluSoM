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

#pragma once

#include "ClosestPoint.h"
#include "PointInsideCheck.h"

namespace SignedDistance
{
	template <typename T>
	double Sphere(T x, T xc, double r)
	{
	    double dis = (x-xc).norm()-r;
	    return dis;
	}

	template <typename T>
	double Cuboid(T xr, T l)
	{
		bool inside = true;
	    double dis = 0.;
	    double maxDis = -1.e32;

	    for (int d=0; d<l.size(); d++)
	    {
	    	double disd = abs(xr(d))-l(d);
	    	if (disd>=0.)
	    	{
	    		dis += disd*disd;
	    		inside = false;
	    	}
	    	else	maxDis = max(maxDis,disd);
	    }
	    double sdis = maxDis;
	    if (!inside) sdis = sqrt(dis);
	    return sdis;
	}

    // under body frame, h is half length
    template <typename T>
	double Cylinder(T xr, double h0, double r0)
	{
		double dis = 0.;
		double h = xr(2);
		double r = sqrt(xr(0)*xr(0) + xr(1)*xr(1));
		double rdif = r-r0;
		double hdif = abs(h)-h0;

		if (hdif<0.)
		{
			if (rdif>0.)	dis = rdif;
			else  			dis = -min(-rdif, -hdif);
		}
		else
		{
			if (rdif>0.)
			{
				xr(2) = 0.;
				dis = sqrt(pow(xr.norm()-r0,2) + hdif*hdif);
			}
			else  			dis = hdif;
		}
		return dis;
	}

	template <typename T0, typename T1>
	double ConvexPolygon2D(T0 x, vector<T0> P, vector<T1> E)
	{
		double dis = numeric_limits<double>::max();
		for (size_t e=0; e<E.size(); ++e)
		{
			T0 A = P[E[e][0]];
			T0 B = P[E[e][1]];
			T0 pe = ClosestPoint::Segment(x,A,B);
			double dise = (pe-x).norm();
			if (dise<dis)	dis = dise;
		}
	    bool xIsInside = PointInsideCheck::ConvexPolygon2D(x, P, E);
	    if (xIsInside) dis *= -1.;
	    return dis;
	}

	template <typename T0, typename T1>
	double Polygon2D(T0 x, vector<T0> P, vector<T1> E)
	{
		double dis = numeric_limits<double>::max();
		for (size_t e=0; e<E.size(); ++e)
		{
			T0 A = P[E[e][0]];
			T0 B = P[E[e][0]];
			T0 pe = ClosestPoint::Segment(x,A,B);
			double dise = (pe-x).norm();
			if (dise<dis)	dis = dise;
		}
	    bool xIsInside = PointInsideCheck::Polygon2D(x, P, E);
	    if (xIsInside) dis *= -1.;
	    return dis;
	}

	template <typename T0, typename T1>
	double Polyhedron3D(T0 x, vector<T0> P, vector<T1> F)
	{
		double dis = numeric_limits<double>::max();
		for (size_t f=0; f<F.size(); ++f)
		{
			T0 cp;
			if (F[f].size()==3)
			{
				cp = ClosestPoint::Triangle(x, P[F[f][0]], P[F[f][1]], P[F[f][2]]);
			}
			else
			{
				cp = ClosestPoint::Polygon3D(x, P, F[f]);
			}
			double disf = (cp-x).norm();
			if (disf<dis)	dis = disf;
		}
		bool xIsInside = PointInsideCheck::Polyhedron(x, P, F);
		if (xIsInside) dis *= -1.;
		return dis;
	}

}