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

#include "../Geometries/GroupedSpheresProperties.h"

class DEM_GROUPED_PARTICLE
{
public:
	DEM_GROUPED_PARTICLE();
	DEM_GROUPED_PARTICLE(size_t tag, vector<DEM_PARTICLE*> lp);
	~DEM_GROUPED_PARTICLE();
	void Conbine();
	void UpdatePosition(double dt);
	void UpdateVelocity(double dt);
	
	DEM_PARTICLE*					P0;
	vector<DEM_PARTICLE*>			Lp;
	vector<Vector3d>				Lr;
	vector<Quaterniond>				Lq;
};

inline DEM_GROUPED_PARTICLE::DEM_GROUPED_PARTICLE(size_t tag, vector<DEM_PARTICLE*> lp)
{
	P0 = new DEM_PARTICLE(tag, Vector3d(0,0,0), 1.);
	Lp = lp;
	Conbine();
}

inline void DEM_GROUPED_PARTICLE::Conbine()
{
	// mass, center
	double m = 0.;
	double vol = 0.;
	Vector3d x (0,0,0);
	for (size_t i=0; i<Lp.size(); ++i)
	{
		x += Lp[i]->X*Lp[i]->Vol;
		m += Lp[i]->M;
		vol += Lp[i]->Vol;
	}
	x /= vol;
	P0->M = m;
	P0->X = x;
	P0->Vol = vol;

	// inertia
	if (isinf(Lp[0]->I(0)))		// for 2d
	{
		double inertia = 0.;
		for (size_t i=0; i<Lp.size(); ++i)
		{
			inertia += Lp[i]->I(2) + Lp[i]->M*pow((Lp[i]->X-x).norm(), 2);
		}
		P0->I = Lp[0]->I;
		P0->I(2) = inertia;

		// double minx = numeric_limits<double>::max();
		// double maxx = -numeric_limits<double>::max();
		// double miny = numeric_limits<double>::max();
		// double maxy = -numeric_limits<double>::max();

		// for (size_t i=0; i<Lp.size(); ++i)
		// {
		// 	minx = min(minx, Lp[i]->X(0)-Lp[i]->R);
		// 	maxx = max(maxx, Lp[i]->X(0)+Lp[i]->R);
		// 	miny = min(miny, Lp[i]->X(1)-Lp[i]->R);
		// 	maxy = max(maxy, Lp[i]->X(1)+Lp[i]->R);
		// }
		// double lx = maxx - minx;
		// double ly = maxy - miny;
	

	}
	else
	{
		cout << "3d combine is not support yet" << endl;
		// 3d need to modify quaternions
		abort();
	}

	// Need to modify Q0, Q, Qf, Qfi for 3d;
	// relative position and orientation under P0
	Matrix3d rotMat = P0->Qfi.toRotationMatrix();
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Vector3d ri = rotMat*(Lp[i]->X - x);
		// cout << "before rotation: " << (Lp[i]->X - x).transpose() << endl;
		// cout << "ri: " << ri.transpose() << endl;
		Lr.push_back(ri);
		Quaterniond qi = Lp[i]->Qfi*P0->Qf;
		Lq.push_back(qi);
	}
	// abort();
}

inline void DEM_GROUPED_PARTICLE::UpdatePosition(double dt)
{
	P0->UpdatePosition(dt);
	Matrix3d rotMat = P0->Qf.toRotationMatrix();
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->X = P0->X + rotMat*Lr[i];
		// rotate inidividual particles
		Lp[i]->Qf = Lq[i]*P0->Qf;
		Lp[i]->Qfi = Lp[i]->Qf.inverse();
	}
}

inline void DEM_GROUPED_PARTICLE::UpdateVelocity(double dt)
{
	Vector3d fc (0,0,0);
	Vector3d tc (0,0,0);
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Vector3d ri = Lp[i]->X - P0->X;
		Vector3d fci = Lp[i]->Fc;
		fc += fci;
		tc += Lp[i]->Tc + ri.cross(fci);
	}
	P0->Fc = fc;
	P0->Tc = tc;

	// Vector3d wb = P0->W;

	P0->UpdateVelocity(dt);
	Vector3d v0 = P0->V;
	Vector3d w0 = P0->Qf._transformVector(P0->W);

	// Vector3d wa = P0->W;
	// if ((wa-wb).norm()>10.)
	// {
	// 	cout << "wb: " << wb.transpose() << endl;
	// 	cout << "wa: " << wa.transpose() << endl;
	// 	cout << "wa-wb: " << (wa-wb).norm() << endl;
	// 	cout << "====================" << endl;
	// }

	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->V = v0 + w0.cross(Lp[i]->X-P0->X);
		Lp[i]->W = Lp[i]->Qfi._transformVector(P0->W);
	}
}