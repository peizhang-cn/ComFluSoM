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

#include "../Mesh/READOBJ.h"
#include "../Geometries/ClosestPoint.h"
#include "../Geometries/PointInsideCheck.h"
#include "../Geometries/SignedDistance.h"
#include "../Geometries/Triangulation.h"
#include "../Geometries/ClosedSurfaceProperties.h"
#include "../Geometries/ConvexityCheck.h"

#ifndef DEM_PARTICLE_H
#define DEM_PARTICLE_H

typedef unordered_map<size_t, Vector3d>	CONVEX_FRICTION_MAP;
typedef unordered_map<size_t, pair<Vector3d, Vector3d>>	CONVEX_META_MAP;
typedef tuple<double, Vector3d, Vector3d>	FRICTION_INFO;
typedef unordered_map<size_t, vector<FRICTION_INFO>>	NONCONVEX_FRICTION_MAP;

class DEM_PARTICLE						    				// class for a single DEM_PARTICLE
{
public:
	DEM_PARTICLE(int tag, const Vector3d& x, double rho);
	// DEM_PATICLE_MOTION
	void UpdatePosition(double dt);
	void UpdateVelocity(double dt);
	void UpdateBox(size_t D);								// Update the range of warpping box and reset cross flags
	void FixV(Vector3d& v);			    					// fix the velocity
	void FixW(Vector3d& w);			    					// fix the angular velocity
	void Fix();												// fix particle with zero velocity
	void UnFix();											// reset all fixs flags
	void SetG(Vector3d& g);			    					// set external acceleration
	void ZeroForceTorque(bool h, bool c);			    	// set zero force and torque
	double CalEnergy();
	Vector3d CalAngularMoment();
	void AddReferencePoint(Vector3d x);
	void MooringBond(size_t ind, double l0, double kn, Vector3d xf);
	Vector3d GetP(size_t i);
	Vector3d Move2BodyFrame(Vector3d x);
	Vector3d Move2GlobalFrame(Vector3d xb);
	// DEM_PARTICLE_ISINSIDE
	bool IsInside(Vector3d x);
	bool IsInsideSphere(Vector3d x);
	bool IsInsideCylinder(Vector3d x);
	bool IsInsideCuboid(Vector3d x);
	bool IsInsidePolygon2D(Vector3d x);
	// DEM_PARTICLE_DISTANCE
	double GetDistance(Vector3d x);
	double GetSignedDistance(Vector3d x);
	double GetSignedDistance2Sphere(Vector3d x);
	double GetSignedDistance2Cuboid(Vector3d x);
	double GetSignedDistance2Cylinder(Vector3d x);
	double GetSignedDistance2Polygon2D(Vector3d x);
	// DEM_PARTICLE_CLOSEST_POINT
	Vector3d GetClosestPoint(Vector3d x);
	Vector3d GetClosestPoint2Sphere(Vector3d x);
	Vector3d GetClosestPoint2Cuboid(Vector3d x);
	Vector3d GetClosestPoint2Cylinder(Vector3d x);
	Vector3d GetClosestPoint2Polygon2D(Vector3d x);
	// DEM_PARTICLE_SET_SHAPE
	void RerangeFaceElementOrder();							// rerange face element order to make sure norm towards to outside, only works for convex shapes
	void SetSphere(double r);								// Change DEM_PARTICLE to Sphere
	void SetDisk2D(double r);								// Change DEM_PARTICLE to 2D disk
	void SetCuboid(double lx, double ly, double lz);		// Change DEM_PARTICLE to Cuboid
	void SetCylinder(double h, double r, Vector3d n);		// Change DEM_PARTICLE to Cylinder
	void SetTriangle2D(vector<Vector3d> ver);
	void SetPolygon2D(vector<Vector3d> ver);
	void SetPolygon3D(vector<Vector3d> ver);
	void SetTetrahedron(vector<Vector3d> ver);
	void SetMetaball2D(double rs, double dis, vector<Vector3d>& metaP, VectorXd& metaK);
	void SetMetaball3D(double rs, double dis, vector<Vector3d>& metaP, VectorXd& metaK);
	void SetFromINP(string fname);

	void InitOBB(size_t d, double e);

    size_t         				ShapeType;                  // ShapeType of DEM_PARTICLE, for 0 is disk2d, for 1 is sphere or 2 is cube etc.
    size_t  					SizeType;					// Size type, 0 is fully solved size, 1 is subgrid size
	size_t         				ID; 				    	// index of DEM_PARTICLE in the list
	int         				Tag;				    	// tag of DEM_PARTICLE
	int         				GroupID;				    // tag of Group
	size_t 						Material;					// material ID is used to find friction coefficient
	size_t 						Nfe;						// Total number of elements in faces
	int							TrajID;						// id for trajectory in the list

	double      				Rho;				    	// density
    double      				R;							// radius of sphere
    double      				Rs;							// radius of sphere factor
	double      				M;					        // mass
	double 						Vol;						// volume
	double      				Kn;					        // normal stiffness (Young's modus if use Hertz model)
	double      				Kt;					        // tangential stiffness
	double						Gn;							// normal viscous damping coefficient (Possion ratio if use Hertz model)
	double						Gt;							// tangential viscous damping coefficient

	Vector3i					Max;						// Max corner of the surrounding box 
	Vector3i 					Min;						// Min corner of the surrounding box 
	Vector3d					BoxL;						// Denmention of the surrounding box
	// Polyhedron
	vector<Vector3d>			P0;				        	// list of point positions at init
	vector<Vector2i>			Edges;				        // list of edges
	vector<VectorXi>			Faces;				        // list of faces
	vector<double>				Parea;						// list of point area for nonconvex metaball
	vector<vector<size_t>>		PNei;						// list of point neighbor for nonconvex metaball
	// Metaball
	vector<Vector3d> 			MetaP0;						// Metaball control points
	VectorXd 					MetaK;						// Metaball coefficients
	double 						MetaC;						// Cricital meta function value for checking cube meta intersection 
	double 						MetaCmin;					// Cricital meta function value for checking cube meta intersection 
	double 						MetaCmax;					// Cricital meta function value for checking cube meta intersection 

	Vector3d					X;				            // position
	Vector3d 					Xt;							// temp position
	Vector3d					Xbr;				        // position before few time step, only used for RWM
	Vector3d					Xb;				            // position before move, only used for DELBM to find refilling LBM nodes
	Vector3d					V;				            // velocity in the center
	Vector3d					W;				            // angular velocity under DEM_PARTICLE frame
	Vector3d					I;				            // inertia under DEM_PARTICLE frame
	Vector3d					G;				        	// Constant body force
	Vector3d					Fh;				        	// Hydro force
	Vector3d					Fc;				        	// Contact force
	Vector3d					Fex;				        // Variable external force that do not need reset to zero
	Vector3d					Th;				        	// Hydro torque under object frame
	Vector3d					Tc;				        	// Contact torque under object frame
	Vector3d					Tex;			        	// Variable external torque under lab frame
	Vector3d					Avb;				        // acceleration of velocity before
	Vector3d					Awb;				        // acceleration of angluar velocity before under object frame
	Vector3d					Xf;							// fixed position
	Vector3d					Vf;				        	// fixed velocity
	Vector3d					Wf;				        	// fixed angylar velocity
	Vector3d					Normal;						// normal direction of wall

    bool        				isRemoved;                	// flag for removed DEM_PARTICLEs
	bool				        isFixV;				    	// flag for fixed translational velocity
	bool        				isFixW;				    	// flag for fixed angular velocity
	bool 						isFixed;					// flag for fixed particle with zero velocity
	bool						isFixTrajectory;			// flag for fixed trajectory
	bool  						isConvex;					// convex or not
	bool  						isBig;						// bigger than bin

    Quaterniond 				Q;				        	// quaternion that describes the rotation
    Quaterniond 				Q0;							// inital quaternion
    Quaterniond					Qf;							// final quaternion (Q0*Q, from object frame to lab frame)
    Quaterniond					Qfi;						// inverse of Qf, used to rotate forces to object frame

    CONVEX_FRICTION_MAP			CFMap;						// convex friction map
    CONVEX_FRICTION_MAP			CFMapt;						// convex friction map updated

    NONCONVEX_FRICTION_MAP		NFMap;						// nonconvex friction map
    NONCONVEX_FRICTION_MAP		NFMapt;						// nonconvex friction map updated

	CONVEX_META_MAP				MMap;						// convex metaball map
	CONVEX_META_MAP				MMapt;						// convex metaball map updated

    Vector3d					OBB[2];						// Oriented Bounding Box

    vector<size_t>				Lb;							// List of bin ID
    vector<pair<size_t, int[3]>>Lc;							// List of particles ID that may in contact
};

#include "DP_ClosestPoint.inl"
#include "DP_Distance.inl"
#include "DP_IsInside.inl"
#include "DP_Motion.inl"
#include "DP_SetShape.inl"
#include "DP_OBB.inl"

inline DEM_PARTICLE::DEM_PARTICLE(int tag, const Vector3d& x, double rho)
{
    ShapeType	= 0;
    SizeType 	= 0;
	ID			= 0;
	Tag			= tag;
	GroupID 	= -1;
	Material 	= 0;
	Nfe 		= 0;
	TrajID		= -1;

	X 		= x;
	Rho		= rho;
	R 		= 0.;
	Rs 		= 0.;
	M 		= 0.;
	Kn	    = 1.0e8;
	Kt		= 1.0e3;
	Gn      = 0.;
	Gt      = 0.;

	Max.setZero();
	Min.setZero();
	BoxL.setZero();

	Xf.setZero();
	V.setZero();
	Vf.setZero();
	W.setZero();
	Wf.setZero();
	G.setZero();
	I << 1., 1., 1.;
	Fh.setZero();
	Fc.setZero();
	Fex.setZero();
	Th.setZero();
	Tc.setZero();
	Tex.setZero();
	Avb.setZero();
	Awb.setZero();
	Normal.setZero();

	Q.w() = 1.;
	Q.vec() << 0.,0.,0.;

	Q0.w() = 1.;
	Q0.vec() << 0.,0.,0.;

	Qf.w() = 1.;
	Qf.vec() << 0.,0.,0.;

	Qfi = Qf.inverse();

	isFixV	= false;
	isFixW	= false;
	isRemoved = false;
	isFixed	= false;
	isFixTrajectory = false;
	isConvex = true;
	isBig 	 = false;

	// CFMap.reserve(30);
	// CFMapt.reserve(30);

	CFMap.clear();
	CFMapt.clear();

	NFMap.clear();
	NFMapt.clear();

	MMap.reserve(10);
	MMapt.reserve(10);
	MMap.clear();
	MMapt.clear();

	// Lb.resize(0);
	Lb.reserve(30);
	// Lc.resize(0);
	Lc.reserve(30);

	P0.resize(0);
	Edges.resize(0);
	Faces.resize(0);
	Parea.resize(0);
	PNei.resize(0);
}

#endif