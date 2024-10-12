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

class MPM_NODE
{
public:
	MPM_NODE();
	MPM_NODE(const Vector3d& x);
	void Reset();
	void NonSlippingBC();
	void SlippingBC(Vector3d& norm);
	void FrictionBC(double dt, Vector3d& n);

	bool 							Actived;					// Flag of actived node

	size_t 							Type;						// 0 for nodes, 1 for gauss quadrature points
	size_t 							ID; 				    	// Index of gird in the list 

	double							M;							// Mass
	double							Vol;						// Volume
	double 							Mu;							// Friction coefficient
	double							SDF;

	double  						JbarDeltaJ;

	Vector3d						X;							// Position
	Vector3d						V;							// Velocity
	Vector3d						Mv;							// Momentum
	Vector3d						F;							// Total force
	Vector3d						Fi;							// Internal force
	Vector3d						Fe;							// External force

	Matrix3d						Stress;						// Stress

	vector<size_t> 					BCTypes;					// Boundary condition type
	vector<Vector3d> 				Norms;						// Normal direction for boundaries
};

inline MPM_NODE::MPM_NODE()
{
	Actived = false;
	Type 	= 0;
	ID 		= 0;
	M 		= 0.;
	Vol 	= 0.;
	SDF 	= 0.;
	X 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	Stress.setZero();
	BCTypes.reserve(3);
	Norms.reserve(3);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline MPM_NODE::MPM_NODE(const Vector3d& x)
{
	Actived = false;
	Type 	= 0;
	ID 		= 0;
	M 		= 0.;
	Vol 	= 0.;
	SDF 	= 0.;
	X 		= x;
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	Stress.setZero();
	BCTypes.reserve(3);
	Norms.reserve(3);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline void MPM_NODE::Reset()
{
	Actived = false;
	M = 0.;
	Vol = 0.;
	JbarDeltaJ = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
	Fi.setZero();
	Fe.setZero();
	Stress.setZero();
}

inline void MPM_NODE::NonSlippingBC()
{
	F.setZero();
	Mv.setZero();
	V.setZero();
}

inline void MPM_NODE::SlippingBC(Vector3d& n)
{	
	// make sure it's compressing the bc then set norm force and momentum to zero
	if (Mv.dot(n)>0.)
	{
		F  = F-F.dot(n)*n;
		Mv = Mv-Mv.dot(n)*n;
		V  = V-V.dot(n)*n;
	}
}

inline void MPM_NODE::FrictionBC(double dt, Vector3d& n)
{
	double mvn = Mv.dot(n);							
	if (mvn>0.)										// 176
	{
		double fn = -mvn/dt;						// 182
		Vector3d mvtV = Mv - mvn*n;
		Vector3d t = mvtV.normalized();				// Tangential vector
		double mvt = mvtV.dot(t);
		double fs = -mvt/dt;						// 184
		double ff = Mu*abs(fn);
		double ft = -min(abs(ff), abs(fs));			// 185

		Vector3d mvb = Mv - F*dt;
		Vector3d mva = (mvt+ft*dt)*t;
		F = (mva - mvb)/dt;
		Mv = mva;
	}
}