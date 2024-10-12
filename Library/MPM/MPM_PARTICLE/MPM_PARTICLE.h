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

struct SHAPE_INFO
{
	size_t		ID	= 0;						// id of node
	double		N	= 0.;						// shape
	Vector3d	Gn	= Vector3d::Zero();			// shape gradient
};

class MPM_PARTICLE
{
public:
	MPM_PARTICLE(size_t n);
	MPM_PARTICLE(size_t n, int tag, const Vector3d& x, double m);
	void SetLinearElastic(double young, double poisson);
	void SetNewtonian(double miu);
	void SetMohrCoulomb(double young, double poisson, double phi, double psi, double c);
	void SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c);

	void UpdatePSize(size_t flag);

	size_t  					CID;

    int 						MID;                        // Type of particle, 0 for elastic 1 for fluid 2 for soil.
	int 						ID; 				    	// Index of particle in the list 
	int 						Tag;				    	// Tag of particle

	double 						M;				            // Mass
	double 						Vol;						// Volume
	double 						Vol0;						// Init Volume

	double						MParas[6];					// Material parameters

	double						JbarDeltaJ;
	double  					JbarDeltaJn;

	double 						P;							// Pressure of fluid or granular materials
// ======================================================================================
	Vector3d					PSize0;						// Vector of half length of particle domain at init
	Vector3d					PSize;						// Vector of half length of particle domain
	Vector3d					Range;						// Vector of influnce range of particle

	Vector3d 					X;				            // Position
	Vector3d 					X0;				            // Init position
	Vector3d					V;							// Velocity
	Vector3d					Vf;							// Fixed velocity
	Vector3d  					B;							// Body force (acceleration)
	Vector3d					Fh;							// Hydro force
	Vector3d					Fc;							// Contact force

	Matrix3d					Stress;						// Stress
	Matrix3d					L;							// Velocity gradient tensor
	Matrix3d					Td;							// Derformation gradient tensor

	bool						FixV;						// Whether the velocity is fixed
	bool						Removed;					// whether this particle is removed

	vector<SHAPE_INFO>			ShapeInfo;					// infomation of shape functions
};

#include "MP_UpdateParticleDomain.inl"
#include "MP_SetLinearElastic.inl"
#include "MP_SetMohrCoulomb.inl"
#include "MP_SetDruckerPrager.inl"
#include "MP_SetNewtonian.inl"

inline MPM_PARTICLE::MPM_PARTICLE(size_t n)
{
    MID		= -1;
	ID		= 0;
	Tag		= 0;
	M 		= 0.;
	X 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fc 		= Vector3d::Zero();

	Stress 	= Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	JbarDeltaJ = 1.;
	JbarDeltaJn = 1.;

	FixV	= false;

	ShapeInfo.reserve(n);
}

inline MPM_PARTICLE::MPM_PARTICLE(size_t n, int tag, const Vector3d& x, double m)
{
    MID		= -1;
	ID		= 0;
	Tag		= tag;
	M 		= m;
	X 		= x;
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fc 		= Vector3d::Zero();

	Range 	= Vector3d::Zero();

	Stress 	= Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	JbarDeltaJ = 1.;
	JbarDeltaJn = 1.;

	FixV	= false;
	Removed	= false;

	ShapeInfo.reserve(n);
}