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

#include "../HEADER.h"
#include "MPM_PARTICLE/MPM_PARTICLE.h"
#include "MPM_NODE/MPM_NODE.h"
#include "MPM_CELL/MPM_CELL.h"
#include "MPM_CELL/MPM_CELL.h"
#include "MPM_SHAPE/MPM_ShapeFunction.h"
#include "../Materials/Materials.h"
#include "../Geometries/SignedDistance.h"

class TOPOGRAPHY
{
	vector<Vector3d>	Vertices;
	vector<VectorXi>	Faces;
};

template<int SType, int D>
class MPM
{
public:
	MPM(): gen(std::random_device()()) {};
	MPM(Vector3d origin, Vector3d lx, double dx, double dt);
	void SetParallel(size_t nproc);
	void FindIndex(size_t n, size_t& i, size_t& j, size_t& k);
	size_t FindIDFromIJK(size_t& i, size_t& j, size_t& k);
	size_t FindIDFromIJK(int& i, int& j, int& k);
	bool CheckIfInDomain(int i, int j, int k);
	void PeriodicNode(int i, int j, int k, int& in, int& jn, int& kn);
	double GetUniformD1();
	/*===================================Functions for init (MPM_init.inl) =====================================================*/		
	void Init(bool useFbar);
	void InitCell();
	/*===================================Functions for init (MPM_AddParticles.inl) =====================================================*/		
	void AddParticle(int tag, Vector3d& x, double m);
	void AddBoxParticles(int tag, Vector3d& x0, Vector3d& x1, double dpx, double rho);
	void DeleteParticles();
	/*===================================Functions for update shape functions (MPM_UpdateShape.inl) =====================================================*/		
	void CalNGN(MPM_PARTICLE* p0);
	void UpdateShapeFunction();
	/*===================================Functions for particle to node (MPM_P2N.inl) =====================================================*/		
	void ParticleToNode();
	/*===================================Functions for particle to cell (MPM_P2C.inl) =====================================================*/		
	void ParticleToCell();
	/*===================================Functions for update on node (MPM_UpdateOnNode.inl) =====================================================*/		
	void CalVOnNode();
	/*===================================Functions for node to particle (MPM_N2P.inl) =====================================================*/		
	void UpdateParticle(MPM_PARTICLE* p0);
	void UpdateStress(MPM_PARTICLE* p0);
	void NodeToParticle();
	/*===================================Functions for double mapping (Reset.inl) =====================================================*/		
	void ResetNodes();
	/*===================================Functions for double mapping (MPM_DoubleMapping.inl) =====================================================*/		
	void NodeToParticleDoubleMapping();
	void CalVOnNodeDoubleMapping();
	/*===================================Functions for boundary conditions (MPM_BC.inl) =====================================================*/		
	void SetNonSlippingBC(size_t n);
	void SetNonSlippingBC(size_t i, size_t j, size_t k);
	void SetSlippingBC(size_t n, Vector3d& norm);
	void SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm);
	void SetFrictionBC(size_t n, double mu, Vector3d& norm);
	void SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm);
	/*===================================Functions for solve (MPM_Solve.inl) =====================================================*/		
	void SolveMUSL(int tt, int ts);
	/*===================================Functions for Input and output (MPM_IO.inl) =====================================================*/		
	void WriteFileH5Particle(string name, int n);
	void WriteFileH5Node(string name, int n);
	void LoadMPMFromH5(string fname, double ratio);

	double 		(*N)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d 	(*GN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	void 		(*NGN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);

	vector <size_t>					LAn;													// List of actived nodes
	vector <size_t>					LAc;													// List of actived nodes
	vector <MPM_PARTICLE*>			Lp;														// List of all MPM particles
	vector <MPM_PARTICLE*>			Lbp;													// List of boundary MPM particles
	vector <MPM_NODE*>				Ln;														// List of all MPM nodes
	vector <MPM_CELL*>				Lcell;													// List of all MPM cells

	vector <Vector3d>				Nei;													// Relative position of neighbor nodes
	vector <Vector3d>				NeiG;													// Shape function gradient of neighbor nodes

	TOPOGRAPHY*						Topology;												// Topology

	bool							Periodic[3];
	bool							ActiveShift;
	bool  							UseFbar;
	bool  							UseCell;
	bool							UseTopology;

	size_t 							Nx;														// Domain size
	size_t 							Ny;
	size_t 							Nz;
	size_t 							Ncz;
	size_t 							Ncy;
	size_t 							Nnode;													// Total number of nodes
	size_t 							Nproc;

	size_t 							ReservedN;												// reserved number of node that influceced by a particle
	size_t  						UpdatePsizeMethod;

	double							Dx;														// Space step
	double 							Dt;														// Time step
	double 							Dc;														// Damping coefficient
	double 							Cs;														// Speed of sound
	double  						Gamma;													// 
	double  						Eta;													// ratio for FLIP and PIC velocities

	Vector3d  						Origin;													// starting point of the mesh
	Vector3d   						Lx;														// size of the mesh
	Vector3d						DomSize;												// size of the mesh + dx (for periodic BC)
	Vector3d  						MinX;													// same as origin
	Vector3d  						MaxX;													// origin + Lx
	Vector3d 						Shift;													// Shift for Galilean invariant
	mt19937 					   	gen;													// Random number generator
};

#include "MPM_BC.inl"
// // #include "MPM_DoubleMapping.inl"
#include "MPM_Init.inl"
#include "MPM_AddParticles.inl"
#include "MPM_IO.inl"
#include "MPM_N2P.inl"
// // #include "MPM_P2C.inl"
#include "MPM_P2N.inl"
#include "MPM_Reset.inl"
#include "MPM_Solve.inl"
#include "MPM_UpdateOnNode.inl"
#include "MPM_UpdateShapeFunction.inl"

template<int SType, int D>
MPM<SType, D>::MPM(Vector3d origin, Vector3d lx, double dx, double dt)
{
	Nproc	= 1;
	Dx 		= dx;
	Dt 		= dt;
	Eta  	= 0.995;
	Dc 		= 0.;
	Cs 		= 0.;
	Gamma	= 7.;

	Origin	= origin;
	Lx 		= lx;

	Nx = ceil(Lx(0)/Dx);
	Ny = ceil(Lx(1)/Dx);
	Nz = ceil(Lx(2)/Dx);

	Periodic[0] = false;
	Periodic[1] = false;
	Periodic[2] = false;

	ActiveShift	= false;
	UseFbar		= false;
	UseCell		= false;
	UseTopology = false;

	Ncz = (Nx+1)*(Ny+1);
	Ncy = (Nx+1);

	Nnode = (Nx+1)*(Ny+1)*(Nz+1);

	Shift.setZero();

	// GIMP
	if (SType==1)
	{
		UpdatePsizeMethod = 0;
		// if (D==1)			NGN =& MPM_ShapeFunction::GIMP1D;
		// else if (D==2) 		NGN =& MPM_ShapeFunction::GIMP2D;
		// else if (D==3)		NGN =& MPM_ShapeFunction::GIMP3D;

		// if (SType==uGIMP)			UpdatePsizeMethod = 0;
		// else if (SType==cpGIMP)		UpdatePsizeMethod = 1;
		// else if (SType==rGIMP)		UpdatePsizeMethod = 2;

		ReservedN = pow(2,D);
		cout << "Using GIMP shape function." << endl;
	}
	else
	{
		cout << "\033[1;31mError: Undefined shape function type! \033[0m\n";		
		exit(0);
	}

	if (D==3)
	{
		Nei  = {   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
	}
	else if (D==2)
	{
		Nei  = {	{ 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
	}

	DomSize = Lx;
	for (size_t d=0; d<D; ++d)		DomSize(d) += Dx;
	MinX = Origin;
	MaxX = DomSize + Origin;
}

// Uniform distribution from -0.5 to 0.5
template<int SType, int D>
inline double MPM<SType, D>::GetUniformD1()
{
	uniform_real_distribution<double>	dis(-0.5,0.5);
	return dis(gen);
}

template<int SType, int D>
inline void MPM<SType, D>::SetParallel(size_t nproc)
{
	Nproc = nproc;
}

// Find index for grid
template<int SType, int D>
inline void MPM<SType, D>::FindIndex(size_t n, size_t& i, size_t& j, size_t& k)
{
	k = n/Ncz;
	j = (n%Ncz)/Ncy;
	i = (n%Ncz)%Ncy;
}

template<int SType, int D>
inline size_t MPM<SType, D>::FindIDFromIJK(size_t& i, size_t& j, size_t& k)
{
	size_t n = i+j*Ncy+k*Ncz;
	return n;	
}

template<int SType, int D>
inline size_t MPM<SType, D>::FindIDFromIJK(int& i, int& j, int& k)
{
	size_t n = i+j*Ncy+k*Ncz;
	return n;
}

template<int SType, int D>
inline bool MPM<SType, D>::CheckIfInDomain(int i, int j, int k)
{
    bool within = false;
    if (i>=0 && i<=(int)Nx && j>=0 && j<=(int)Ny && k>=0 && k<=(int)Nz) within = true;
    return within;
}

template<int SType, int D>
inline void MPM<SType, D>::PeriodicNode(int i, int j, int k, int& in, int& jn, int& kn)
{
    in = i; jn = j; kn = k;
    if (Periodic[0])    in = (i+Nx+1)%(Nx+1);
    else                in = i;
    if (Periodic[1])    jn = (j+Ny+1)%(Ny+1);
    else                jn = j;
    if (Periodic[2])    kn = (k+Nz+1)%(Nz+1);
    else                kn = k;
}