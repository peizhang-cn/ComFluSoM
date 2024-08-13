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
 
// This part of code is used to calculate mass center, volume and inertia tensor of polyhedron

#pragma once

#ifndef PARTICLE_PROPERTIES_H
#define PARTICLE_PROPERTIES_H

// "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates"
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
// only works for triangle mesh
void CalTriangleMeshProperties(vector<VectorXi>& F, vector<Vector3d>& P, double rho, double& vol, double& m, Vector3d& xc, Matrix3d& inertia)
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

#endif