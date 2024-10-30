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

namespace ClosedSurfaceProperties
{
// calculate the area of a polygon in 2d
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

	void CalTetrahedronInertia(Vector3d& P1, Vector3d& P2, Vector3d& P3, Vector3d& P4, Vector3d& xc, double& m, Matrix3d& inertia)
	{
		Vector3d p1 = P1-xc;
		Vector3d p2 = P2-xc;
		Vector3d p3 = P3-xc;
		Vector3d p4 = P4-xc;

		Vector3d a2;
		for (size_t i=0; i<3; ++i)
		{
			a2(i) = p1(i)*p1(i)+p1(i)*p2(i)+p2(i)*p2(i)+p1(i)*p3(i)+p2(i)*p3(i)+p3(i)*p3(i)+p1(i)*p4(i)+p2(i)*p4(i)+p3(i)*p4(i)+p4(i)*p4(i);
		}
		double a = 0.1*m*(a2(1)+a2(2));
		double b = 0.1*m*(a2(0)+a2(2));
		double c = 0.1*m*(a2(0)+a2(1));

		double aa = 2.*p1(1)*p1(2)+p2(1)*p1(2)+p3(1)*p1(2)+p4(1)*p1(2)+p1(1)*p2(2)+2.*p2(1)*p2(2)+p3(1)*p2(2)+p4(1)*p2(2)+p1(1)*p3(2)+p2(1)*p3(2)+2.*p3(1)*p3(2)+p4(1)*p3(2)+p1(1)*p4(2)+p2(1)*p4(2)+p3(1)*p4(2)+2.*p4(1)*p4(2);
		double bb = 2.*p1(0)*p1(2)+p2(0)*p1(2)+p3(0)*p1(2)+p4(0)*p1(2)+p1(0)*p2(2)+2.*p2(0)*p2(2)+p3(0)*p2(2)+p4(0)*p2(2)+p1(0)*p3(2)+p2(0)*p3(2)+2.*p3(0)*p3(2)+p4(0)*p3(2)+p1(0)*p4(2)+p2(0)*p4(2)+p3(0)*p4(2)+2.*p4(0)*p4(2);
		double cc = 2.*p1(0)*p1(1)+p2(0)*p1(1)+p3(0)*p1(1)+p4(0)*p1(1)+p1(0)*p2(1)+2.*p2(0)*p2(1)+p3(0)*p2(1)+p4(0)*p2(1)+p1(0)*p3(1)+p2(0)*p3(1)+2.*p3(0)*p3(1)+p4(0)*p3(1)+p1(0)*p4(1)+p2(0)*p4(1)+p3(0)*p4(1)+2.*p4(0)*p4(1);

		aa *= 0.05*m;
		bb *= 0.05*m;
		cc *= 0.05*m;

		inertia(0,0) = a;
		inertia(1,1) = b;
		inertia(2,2) = c;
		inertia(0,1) = inertia(1,0) = -cc;
		inertia(0,2) = inertia(2,0) = -bb;
		inertia(1,2) = inertia(2,1) = -aa;
	}

	void CalTriangleMeshProperties(vector<Vector3d>& P, vector<VectorXi>& F, double rho, double& vol, double& m, Vector3d& xc, Matrix3d& inertia)
	{
		vol = 0.;
		m = 0.;
		xc.setZero();
		inertia.setZero();

		Vector3d pr (0.,0.,0.);				// reference point
		for (size_t i=0; i<P.size(); ++i)
		{
			pr += P[i];
		}
		pr /= P.size();

		for (size_t i=0; i<F.size(); ++i)
		{
			Vector3d p0 = P[F[i](0)];	// vertice of the face
			Vector3d p1 = P[F[i](1)];
			Vector3d p2 = P[F[i](2)];
			Vector3d pc = (p0+p1+p2)/3.;	// triangle center
			Vector3d norm = (p1-p0).cross(p2-p0);	// outside normal

			double voli = abs((p0-pr).dot((p1-pr).cross(p2-pr)))/6.;	// volume
			Vector3d xci = (p0+p1+p2+pr)/4.;							// center
			
			double sign = 1.;
			if ((pc-pr).dot(norm)<0.)	sign = -1.;
			voli *= sign;					// signed volume

			double mi = rho*abs(voli);
			Matrix3d inertiai;
			CalTetrahedronInertia(p0, p1, p2, pr, pr, mi, inertiai);
			vol += voli;
			xc += xci*voli;
			inertia += sign*inertiai;
		}
		xc /= vol;
		m = rho*vol;
		// move inertia to mass center
		Vector3d xcr = xc-pr;
		inertia(0,0) -= m*(xcr(1)*xcr(1) + xcr(2)*xcr(2));
		inertia(1,1) -= m*(xcr(0)*xcr(0) + xcr(2)*xcr(2));
		inertia(2,2) -= m*(xcr(0)*xcr(0) + xcr(1)*xcr(1));
		inertia(0,1) += m*xcr(0)*xcr(1);	inertia(1,0) = inertia(0,1);
		inertia(0,2) += m*xcr(0)*xcr(2);	inertia(2,0) = inertia(0,2);
		inertia(1,2) += m*xcr(1)*xcr(2);	inertia(2,1) = inertia(1,2);
	}

	void CalPolygon2DProperties(const vector<Vector3d>& ver, double rho, double& area, double& m, Vector3d& xc, double& Iz)
	{
		size_t nv = ver.size();
		double cx = 0., cy = 0., jx = 0., jy = 0., sa = 0.;
		Vector3d ref = ver[0];  // Reference point for improved numerical stability

		for (size_t i = 0; i < nv; ++i)
		{
			size_t j = (i+1) % nv;
			Vector3d vi = ver[i] - ref, vj = ver[j] - ref;
			double fact = (vi(0)*vj(1) - vj(0)*vi(1));
			cx += (vi(0)+vj(0))*fact;
			cy += (vi(1)+vj(1))*fact;
			jx += (vi(1)*vi(1) + vi(1)*vj(1) + vj(1)*vj(1))*fact;
			jy += (vi(0)*vi(0) + vi(0)*vj(0) + vj(0)*vj(0))*fact;
			sa += fact;
		}

		area = abs(0.5*sa);
		m = rho * area;

		if (sa != 0) {  // Avoid division by zero
			xc(0) = cx/(3.*sa) + ref(0);
			xc(1) = cy/(3.*sa) + ref(1);
			xc(2) = 0.;
			Iz = rho*(jx+jy)/12. - m*((xc-ref).squaredNorm());
		} else {
			xc = Vector3d::Zero();
			Iz = 0;
			cout << "Warning: Area is zero!" << endl;
			abort();
		}
	}
}