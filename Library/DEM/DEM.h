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

#include "../HEADER.h"
#include "DEM_PARTICLE/DEM_PARTICLE.h"
#include "../LinkedCell/BINS_LC.h"
#include "DEM_MAPS.h"
// #include "GJK.h"
#include "TRAJECTORY.h"
#include "META_SURFACE_MESH.h"
#include "METABALL.h"
#include "PARTICLE_PROPERTIES.h"
#include "PARTICLE_PROPERTIES_CONVEX.h"

#ifndef DEM_H
#define DEM_H

struct MOORING_BOND
{
	size_t I;
	size_t P;
	double Kn;
	double L;
	Vector3d Xf;
};

struct CONTACT_INFO
{
	double 		Delta = 0.;					// overlaps
	Vector3d	Xc = Vector3d::Zero();		// contact point
	Vector3d	Nc = Vector3d::Zero();		// contact normal for first particle
};

struct CONTACT_PARA
{
	double Kn = 0.;
	double Gn = 0.;
	double Kt = 0.;
	double Gt = 0.;
};

class DEM
{
public:
	DEM(Vector3d origin, Vector3d l, string cmtype, double dt);
	~DEM();
	void SetPeriodic(bool bx, bool by, bool bz);															// Set periodic bounary conditions
	void InitBinSystem(VectorXi n);																			// Generate bin system
	void SetParallel(size_t nproc);																			// Set number of processor to use
	void SetMaterialCoef(size_t i, size_t j, double cr, double mus, double mud);
	void ActiveMetaball();																					// Active metaball functions
	void ActiveREV();																						// Active REV
	void SetLubrication(double hn, double viscosity);
	void SetG(Vector3d& g);																					// Set gravity
	void SetTrajectory(DEM_PARTICLE* pi, vector<Vector3d>& x, vector<Vector3d>& v, vector<Vector3d>& w);	// Set trajectory
	void SetTrajectory(DEM_PARTICLE* pi, size_t tid);														// Set trajectory
	void DeleteParticles();
	void MappingREV(bool l2e);
	void SetMooringBond(size_t p, size_t ind, double l0, double kn, Vector3d xf);
	void ApplyMooringBond();
	size_t Key(int i, int j);
	/*===================================Functions for add particles (DEM_ADD_PARTICLES.h) =====================================================*/
	void AddDisk2D(int tag, double r, Vector3d& x, double rho);								// Add a 2d disk
	void AddSphere(int tag, double r, Vector3d& x, double rho);								// Add a sphere
	void AddMetaball(int tag, double dis, vector<Vector3d>& metaP, VectorXd& metaK, double rho);
	void AddCuboid(int tag, Vector3d& l, Vector3d& x, double rho);
	void AddCylinder(int tag, double h, double r, Vector3d& n, Vector3d& x, double rho);	// Add a cylinder
	// void AddTetrahedron(int tag, vector<Vector3d> ver, double rho);
	// void AddTriangle2D(int tag, vector<Vector3d> ver, double rho);
	// void AddPolygon2D(int tag, vector<Vector3d> ver, double rho);
	// void AddPolygon3D(int tag, double r, vector<Vector3d> ver, double rho);
	void AddNSpheres(int tag, size_t seed, size_t np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho);	// Add np spheres
	void AddNMetaballs(int tag, size_t seed, size_t np, Vector3d& x0, Vector3d& x1, double rho, double dis, vector<Vector3d>& metaP, VectorXd& metaK);
	/*===================================Functions for contact detection (DEM_CONTACT_DETECTION.h) =====================================================*/		
	void LinkedCell(bool firststep);
	/*===================================Functions for calculate Contact parameters (DEM_CONTACT_PARAMETERS.h) =====================================================*/	
	double EffectiveValue(double ai, double aj);													// Calculate effective values for contact force
	void LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, CONTACT_PARA& para);					// Calculate parameters for linear contact model
	void LinearContactParaFromCr(DEM_PARTICLE* pi, DEM_PARTICLE* pj, CONTACT_PARA& para);
	void HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, CONTACT_PARA& para);	// Calculate parameters for Hertz contact model
	void ContactPara(int cmType, DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, CONTACT_PARA& para);
	/*===================================Functions for calculate Contact between sphere and cylinder (DEM_CONTACT_SPHERE_CYLINDER.h) =====================================================*/	
	void Sphere2Cylinder(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, vector<CONTACT_INFO>& lci);
	/*===================================Functions for calculate Contact between sphere and cylinder (DEM_CONTACT_SPHERE_CUBIOD.h) =====================================================*/	
	void Sphere2Cuboid(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, vector<CONTACT_INFO>& lci);
	/*===================================Functions for calculate Contact between convex Polyhedras (DEM_CONTACT_CONVEX_POLYHEDRA.h) =====================================================*/	
	void Polyhedra2PolyhedraConvex(Vector3d xi, Vector3d xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, vector<CONTACT_INFO>& lci);
	void Sphere2PolyhedraConvex(Vector3d xi, Vector3d xj, DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, vector<CONTACT_INFO>& lci);
	/*===================================Functions for calculate Contact between convex metaballs (DEM_CONTACT_CONVEX_METABALL.h) =====================================================*/		
	void Metaball2Metaball(Vector3d xi, Vector3d xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, DEM_MAPS* maps, Vector3d& finalP, double& delta, Vector3d& n, Vector3d& cp);
	void Metaball2Wall(DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, DEM_MAPS* maps, Vector3d& finalP, double& delta, Vector3d& n, Vector3d& cp);
	/*===================================Functions for contact (DEM_CONTACT.h) =====================================================*/		
	void Contact2P(int cmType, DEM_PARTICLE* pi, DEM_PARTICLE* pj, bool forP, int peri[3]);
	void PeriodiParticlePosition2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, bool forP, int peri[3], Vector3d& Xi, Vector3d& Xj);
	void ContactInformation2P(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, vector<CONTACT_INFO>& lci);
	void FrictionTangentForce(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, double kt, double gt, Vector3d& cp, Vector3d& n, Vector3d& fn, Vector3d& ft);
	// void RollingResistance(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kr, double gr, Vector3d& n, Vector3d& fn, Vector3d& xir, Vector3d& armr);
	void ContactForce2P(int cmType, Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, vector<CONTACT_INFO>& lci);
	void Contact(bool writeFc, int n, bool show);
	/*===================================Functions for contact (DEM_MOVE.h) =====================================================*/		
	void UpdatePositionGlobal(size_t t, double dt, bool resetHydro);						// Update the particle position by velocity verlet
	void UpdateVelocityGlobal(double dt);													// Update the particle velocity by velocity verlet
	void ZeroForceTorqueGlobal(bool h, bool c);												// Reset forces to zero
	/*===================================Functions for solve (DEM_SOLVE.h) =====================================================*/		
	void SolveOneStep(size_t t, size_t ts, vector<double>& times);
	void Solve(int tt, int ts, bool writefile);
	/*===================================Functions for Input and output (DEM_IO.h) =====================================================*/		
	void WriteFileH5(int n);
	void WriteContactForceFileH5(int n);

    string 							CMType;													// Contact force model type

    bool 							Periodic[3];											// Flag for periodic conditions

    bool 							BinExist;												// True if bin is set
    bool  							UsingMetaball;											// True if using metaball
    bool  							UsingREV;

    size_t 							D;														// Dimension
    size_t 							Np;														// Total number of points in the domain
    size_t 							Nf;														// Total number of faces in the domain
    size_t 							Nproc;													// Number of processors
    size_t  						MetaLevel;												// Level of cell that used to find inital point

    int  							CMTypeIndex;

    double							Lx;														// Domain size
    double							Ly;														// Domain size
    double							Lz;														// Domain size	
    double 							Dt;														// Time step

    Vector3d  						Origin;
    Vector3d  						DomMax;
    Vector3d  						DomSize;

	double 							Hn;														// Lubrication cut off distance for fluid
	double 							Viscosity;												// Dynamic viscosity fluid

	bool  							MaterialTable[10][10];									// bool to indicate if materials are defined
	double 							CrTable[10][10];										// Coefficient of restitution table

	double 							FsTable[10][10];										// Friction coefficient (static) table
	double 							FdTable[10][10];										// Friction coefficient (dynamic) table
	double 							RsTable[10][10];										// Rolling friction coefficient (static) table
	double 							RdTable[10][10];										// Rolling Friction coefficient (dynamic) table

	vector<size_t>  				Lmaterial;

	vector<Vector3d> 				UnitSpherePoints;										// Unit sphere mesh points for metaball
	vector<VectorXi> 				UnitSphereFaces;										// Unit sphere mesh faces for metaball

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < DEM_PARTICLE* >		Lg;														// List of groups
	vector<size_t>					Lw;														// List of particle that are larger than bin, typical a wall

	vector<vector<size_t>> 			Lfid;													// List of global face id to particle id
	vector<TRAJECTORY*>				Ltraj;													// List of particle trajectories

	vector<MOORING_BOND*>			Lmb;

	double  						VolRev0;												// Initial REV volume
	double  						VolRev;													// REV volume
	Matrix3d						Td;														// Mapping from original config to current
	Matrix3d						Tdi;													// Mapping from current to original config
	Matrix3d  						Stress;

	BINS_LC*						Bins;													// Bin system for Linked Cell method
};

#include "DEM_AddParticle.inl"
#include "DEM_Contact.inl"
#include "DEM_ContactDetection.inl"
#include "DEM_ContactPara.inl"
#include "DEM_ConvexMetaball.inl"
#include "DEM_ConvexPolyhedra.inl"
#include "DEM_IO.inl"
#include "DEM_Move.inl"
#include "DEM_Solve.inl"
#include "DEM_SphereCubiod.inl"
#include "DEM_SphereCylinder.inl"

inline DEM::DEM(Vector3d origin, Vector3d l, string cmtype, double dt)
{
	BinExist = false;
	UsingMetaball = false;

	UsingREV = false;
	Td.setIdentity();
	Tdi.setIdentity();

	Origin = origin;
	DomMax = origin+l;
	DomSize = l;
	Lx = DomSize(0);	Ly = DomSize(1);	Lz = DomSize(2);

	Dt = dt;
	Np = 0;
	Nf = 0;

	D = (Lz==0.) ? 2:3;

	Nproc = 1;
	MetaLevel = 4;

	Periodic[0] = true;	Periodic[1] = true;	Periodic[2] = true;

	Hn = 0.;
	Viscosity = 0.;

	CMType = cmtype;

	if (CMType=="LINEAR")
	{
		CMTypeIndex = 0;
		cout << "Using Linear contact model." << endl;
	}
	else if (CMType=="LINEAR_CR")
	{
		CMTypeIndex = 1;
		cout << "Using Linear contact model and damping with Cr." << endl;
	}
	else if (CMType=="HERTZ")
	{
		CMTypeIndex = 2;
		cout << "Using Hertz Contact model." << endl;
	}
	else
	{
		cout << "Undefined contact model!" << endl;
		abort();
	}
	Lmaterial.resize(0);
	Lp.resize(0);
	Lg.resize(0);
	Lw.resize(0);
	Lfid.resize(0);
	Ltraj.resize(0);
	Lmb.resize(0);
	// Set friction coefficient to zero
	for (size_t i=0; i<10; ++i)
	for (size_t j=0; j<10; ++j)
	{
		MaterialTable[i][j] = false;
		CrTable[i][j] = 1.;
		FsTable[i][j] = 0.;
		FdTable[i][j] = 0.;
		RsTable[i][j] = 0.;
		RdTable[i][j] = 0.;
	}
}

inline void DEM::SetPeriodic(bool bx, bool by, bool bz)
{
	if (BinExist)
	{
		cout << "\033[1;31mError: Need to set periodic boundaries before generate bin system! \033[0m\n";		
		exit(0);			
	}
	if (!bx) Periodic[0] = false;
	if (!by) Periodic[1] = false;
	if (!bz) Periodic[2] = false;
}

inline void DEM::SetParallel(size_t nproc)
{
	Nproc = nproc;
	if (BinExist)	Bins->SetParallel(nproc);
	else
	{
		cout << "SetParallel" << endl;
		cout << "Need to have bin first!" << endl;
		abort();
	}
}

inline void DEM::SetMaterialCoef(size_t i, size_t j, double cr, double mus, double mud)
{
	if (!MaterialTable[i][j])
	{
		CrTable[i][j] = CrTable[j][i] = cr;
		FsTable[i][j] = FsTable[j][i] = mus;
		FdTable[i][j] = FdTable[j][i] = mud;
		MaterialTable[i][j] = MaterialTable[j][i] = true;
	}
	else
	{
		cout << "\033[33m" << "Noting: Material Coef already defined, i: " << i << " and material j: " << j << endl;
		cout << "\033[" << 0 << "m";
		// abort();
	}
}

inline void DEM::InitBinSystem(VectorXi n)
{
	Vector3d l (Lx, Ly, Lz);
	Bins = new BINS_LC(Origin, l, n, Periodic[0], Periodic[1], Periodic[2]);
	BinExist = true;
}

// inline void DEM::ActiveMetaball()
// {
// 	vector<Vector3d> fn;
// 	ReadOBJ("UnitSphere.obj", UnitSpherePoints, UnitSphereFaces, fn);
// 	cout << "UnitSpherePoints: " << UnitSpherePoints.size() << endl;
// 	cout << "Metaball actived." << endl;
// 	UsingMetaball = true;
// }

// inline void DEM::ActiveREV()
// {
// 	UsingREV = true;
// 	if (D==2)		VolRev0 = Lx*Ly;
// 	else if (D==3)	VolRev0 = Lx*Ly*Lz;
// 	VolRev = VolRev0;
// }

// Map a pair of integer to one integer key for hashing
// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
inline size_t DEM::Key(int i, int j)
{
	return (i+j+1)*(i+j)/2+j;
}

// inline void DEM::MappingREV(bool l2e)
// {
// 	if (l2e)
// 	{
// 		#pragma omp parallel for schedule(static) num_threads(Nproc)
// 		for (size_t i=0; i<Lp.size(); ++i)
// 		{
// 			Lp[i]->Xt = Lp[i]->X;
// 			Lp[i]->X = Tdi*Lp[i]->Xt;
// 		}
// 	}
// 	else
// 	{
// 		#pragma omp parallel for schedule(static) num_threads(Nproc)
// 		for (size_t i=0; i<Lp.size(); ++i)
// 		{
// 			Lp[i]->X = Lp[i]->Xt;
// 		}		
// 	}
// }

inline void DEM::SetMooringBond(size_t i, size_t p, double l0, double kn, Vector3d xf)
{
	MOORING_BOND* mb = new MOORING_BOND();
	mb->I = i;
	mb->P = p;
	mb->L = l0;
	mb->Kn = kn;
	mb->Xf = xf;
	Lmb.push_back(mb);
}

inline void DEM::ApplyMooringBond()
{
	for (size_t m=0; m<Lmb.size(); ++m)
	{
		size_t i = Lmb[m]->I;
		size_t p = Lmb[m]->P;
		Lp[i]->MooringBond(p, Lmb[m]->L, Lmb[m]->Kn, Lmb[m]->Xf);
	}
}

inline void DEM::SetG(Vector3d& g)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->SetG(g);
	}
}

inline void DEM::SetTrajectory(DEM_PARTICLE* pi, vector<Vector3d>& x, vector<Vector3d>& v, vector<Vector3d>& w)
{
	Ltraj.push_back(new TRAJECTORY(x,v,w));
	size_t tid = Ltraj.size()-1;
	pi->TrajID = tid;
	pi->isFixTrajectory = true;
	pi->isFixV = true;
	pi->isFixW = true;
	if ((pi->X-Ltraj[tid]->X[0]).norm()>1.e-20)
	{
		cout << "the starting point of the particle is not equal to the starting point of this Trajectory" << endl;
		abort();
	}
}

inline void DEM::SetTrajectory(DEM_PARTICLE* pi, size_t tid)
{
	if (tid>=Ltraj.size())
	{
		cout << "the trajectory does not exist!" << endl;
		abort();
	}
	pi->TrajID = tid;
	pi->isFixTrajectory = true;
	pi->isFixV = true;
	pi->isFixW = true;
	if ((pi->X-Ltraj[tid]->X[0]).norm()>1.e-20)
	{
		cout << "the starting point of the particle is not equal to the starting point of this Trajectory" << endl;
		abort();
	}	
}

inline void DEM::SetLubrication(double hn, double viscosity)
{
	Hn = hn;
	Viscosity = viscosity;
}

inline void DEM::ContactPara(int cmType, DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, CONTACT_PARA& para)
{
	switch(cmType)
	{
	case 0:
		LinearContactPara(pi, pj, para);
		break;
	case 1:
		LinearContactParaFromCr(pi, pj, para);
		break;
	case 2:
		HertzContactPara(pi, pj, delta, para);
		break;
	}
}

// inline void DEM::DeleteParticles()
// {
// 	vector <DEM_PARTICLE*>	Lpt;
// 	Lpt.resize(0);

// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		if (!Lp[p]->removed)	Lpt.push_back(Lp[p]);
// 	}
// 	Lp = Lpt;

// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		Lp[p]->ID = p;
// 	}
// }

#endif