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

#ifndef __CLOESTPOINT_H__
#define __CLOESTPOINT_H__

namespace ClosestPoint
{
	// find closest point on a sphere to point x
	template <typename T>
	T SphereClosestPoint(T x, T xc, double r)
	{
		T vec = x-xc;
		vec.normalize();
		T cp = x + r*vec;
	    return cp;
	}

	// find closest point on a sphere to point x
	template <typename T>
	T CylinderClosestPoint(T xr, double h0, double r0)
	{
		double h = xr(2);
		double r = sqrt(xr(0)*xr(0) + xr(1)*xr(1));

	    T cpr;
	    cpr.setZero();

	    if (r>r0)   cpr = r0/r*xr;
	    else        cpr = xr;

	    if (h>h0)           cpr(2) = h0;
	    else if (h<-h0)     cpr(2) = -h0;
	    else                cpr(2) = xr(2);
		return cpr;
	}

	// find closest point on a segment to point x (works for 2D and 3D, Vector2d, T)
	template <typename T>
	T SegmentClosestPoint(T& x, T& A, T& B)
	{
		T AB = B-A;
		double k = (x-A).dot(AB)/AB.squaredNorm();
		T p = A;
		if (k>1.)		p = B;
		else if (k>0.)	p += k*AB;
		return p;
	}

	// find closest point on a Triangle to point x, only for 3D (also works for 2D but not optimized for efficiency)
	template <typename T>
	T TriangleClosestPoint(T& x, T& A, T& B, T& C)
	{
		T AC = C-A;
		T AB = B-A;
		T XA = A-x;
		T ABC = AB.cross(AC);
		T p;

		if (XA.dot(ABC.cross(AC))<0.)		p = SegmentClosestPoint(x, A, C);
		else if (XA.dot(AB.cross(ABC))<0.)	p = SegmentClosestPoint(x, A, B);
		else
		{
			T XB = B-x;
			T BC = C-B;
			if (XB.dot(BC.cross(ABC))<0.)	p = SegmentClosestPoint(x, B, C);
			else  							p = x+ABC.dot(XA)/ABC.squaredNorm()*ABC;
		}
		return p;
	}

	// find closest point on a Quadrilateral to point x, only for 3D
	template <typename T>
	T QuadrilateralClosestPoint(T& x, T& A, T& B, T& C, T& D)
	{
		T AD = D-A;
		T AB = B-A;
		T XA = A-x;
		T ABD = AB.cross(AD);
		T p;

		if (XA.dot(ABD.cross(AD))<0.)		p = SegmentClosestPoint(x, A, D);
		else if (XA.dot(AB.cross(ABD))<0.)	p = SegmentClosestPoint(x, A, B);
		else
		{
			T CB = B-C;
			T CD = D-C;
			T XC = C-x;
			if (XC.dot(ABD.cross(CB))<0.)		p = SegmentClosestPoint(x, B, C);
			else if (XC.dot(CD.cross(ABD))<0.)	p = SegmentClosestPoint(x, C, D);
			else								p = x+ABD.dot(XA)/ABD.squaredNorm()*ABD;
		}
		return p;
	}

	template <typename T>
	T Polygon2DClosestPoint(T& x, vector<T> P, vector<Vector2i> E)
	{
		double dis = 1.e30;
		T p;
		for (size_t e=0; e<E.size(); ++e)
		{
			T A = P[E[e](0)];
			T B = P[E[e](1)];
			T pe = SegmentClosestPoint(x,A,B);
			double dise = (pe-x).norm();
			if (dise<dis)
			{
				dis = dise;
				p = pe;
			}
		}
		return p;
	}

	template <typename T>
	T Polygon3DClosestPoint(T& x, T& norm, vector<T>& P, vector<Vector2i>& E)
	{
		double dis = 1.e30;
		T p;
		T n = norm;			// normal

		bool isFace = true;
		for (size_t e=0; e<E.size(); ++e)
		{
			T A = P[E[e](0)];
			T B = P[E[e](1)];

			T ne = (B-A).cross(n);
			double sign = (x-A).dot(ne);
			if (sign>0.)
			{
				isFace = false;
				T pe = SegmentClosestPoint(x,A,B);
				double dise = (pe-x).norm();
				if (dise<dis)
				{
					dis = dise;
					p = pe;
				}
			}
		}
		if (isFace)
		{
			// cout << "isface" << endl;
			p = x-(x-P[E[0](0)]).dot(n)*n;
		}

		// p = x-(x-P[E[0](0)]).dot(n)*n;

		return p;
	}

	template <typename T>
	T Polygon3DClosestPoint(T& x, vector<T>& P)
	{
		double dis = 1.e30;
		T p;
		T n = (P[1]-P[0]).cross(P[2]-P[1]);			// normal
		n.normalize();

		size_t np = p.size();
		bool isFace = true;
		for (size_t i=0; i<np; ++i)
		{
			T A = P[i];
			T B = P[(i+1)%np];

			T ne = (B-A).cross(n);
			double sign = (x-A).dot(ne);
			if (sign>0.)
			{
				isFace = false;
				T pe = SegmentClosestPoint(x,A,B);
				double dise = (pe-x).norm();
				if (dise<dis)
				{
					dis = dise;
					p = pe;
				}
			}
		}
		if (isFace)
		{
			// cout << "isface" << endl;
			p = x-(x-P[0]).dot(n)*n;
		}
		return p;
	}

	template <typename T>
	T Polygon3DClosestPoint(T& x, vector<T>& P, VectorXi face)
	{
		T p;

		size_t ne = face.size();
		vector<Vector2i> edges(0);
		for (size_t i=0; i<ne; ++i)
		{
			Vector2i edgei;
			edgei(0) = face(i);
			edgei(1) = face((i+1)%ne);
			edges.push_back(edgei);
		}
		T n = (P[edges[0](1)]-P[edges[0](0)]).cross(P[edges[1](1)]-P[edges[1](0)]);
		n.normalize();
		p = Polygon3DClosestPoint(x, n, P, edges);
		return p;
	}

	template <typename T>
	T CuboidClosestPoint(T xr, T l)
	{
		bool inside = true;
		T cp = xr;
		// cp.setZero();

	    for (int d=0; d<l.size(); d++)
	    {
	    	if (xr(d)-l(d)>0.)
	    	{
	    		cp(d) = l(d);
	    		inside = false;
	    	}
	    	else if (xr(d)+l(d)<0.)
	    	{
	    		cp(d) = -l(d);
	    		inside = false;
	    	}
	    }

		double minDis = 1.e32;
		int side = 0;
		int ori = 0;

	    if (inside)	
	    {
	    	for (int d=0; d<l.size(); d++)
	    	{
	    		double dis = l(d)-abs(xr(d));
	    		if (dis<minDis)
	    		{
	    			side = (xr(d)>0.) ? 1:-1;
	    			minDis = dis;
	    			ori = d;
	    		}
	    	}
	    	cp = xr;
	    	cp(ori) = side*l(ori);
	    }

		return cp;
	}
}

#endif