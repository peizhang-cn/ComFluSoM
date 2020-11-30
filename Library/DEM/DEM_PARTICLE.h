/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#include <PARTICLE_PROPERTIES.h>
#include <METABALL.h>

class DEM_PARTICLE						    				// class for a single DEM_PARTICLE
{
public:
	DEM_PARTICLE(int tag, const Vector3d& x, double rho);
	void Set(double kn);			    					// set physical parameters of the DEM_PARTICLEs
	void SetG(Vector3d& g);			    					// set external acceleration
	void VelocityVerlet(double dt, bool movePoints);		// move the DEM_PARTICLE based on Velocity Verlet intergrator
	void UpdatePosition(double dt, bool movePoints, bool resetHydro);
	void UpdateVelocity(double dt);
	void UpdateBox(size_t D);								// Update the range of warpping box and reset cross flags
	void Constrain(int d, double dt);						// Constrain on d direction with v velocity, applied after VelocityVerlet
	void FixV(Vector3d& v);			    					// fix the velocity
	void FixW(Vector3d& w);			    					// fix the angular velocity
	void Fix();												// fix particle with zero velocity
	void UnFix();											// reset all fixs flags
	double CalEnergy();
	Vector3d CalAngularMoment();
	void ZeroForceTorque(bool h, bool c);			    	// set zero force and torque
	void RerangeFaceElementOrder();							// rerange face element order to make sure norm towards to outside, only works for convex shapes
	void SetPlane(Vector3d& n);								// Change DEM_PARTICLE to infinite plane
	void SetDrum(double rd, Vector3d& xd);
	void SetSphere(double r);								// Change DEM_PARTICLE to Sphere
	void SetDisk2D(double r);								// Change DEM_PARTICLE to 2D disk
	void Set2DPolynomialParticle(VectorXd& coef);			// Change DEM_PARTICLE to 2D Polynomial Particle
	void SetCuboid(double lx, double ly, double lz);		// Change DEM_PARTICLE to Cuboid
	void SetTetrahedron(vector<Vector3d> ver);
	void UpdateCoef();
	void SetMetaball(double rs, vector<Vector3d>& metaP, VectorXd& metaK, vector<Vector3d>& unitSphereP, vector<VectorXi>& unitSphereF);
	// void DistanceToSurface(Vector3d& x);

	size_t 						PID;						// processor id
	int							CID;						// bin id for linked cell
	int							CIDt;						// tempt bin id for linked cell

    int         				Type;                       // Type of DEM_PARTICLE, for 0 is disk2d, for 1 is sphere or 2 is cube etc.
	int         				ID; 				    	// index of DEM_PARTICLE in the list 
	int         				Tag;				    	// tag of DEM_PARTICLE
	int         				Group;				    	// tag of Group
	int 						MID;						// material ID is used to find friction coefficient
	int 						Nfe;						// Total number of elements in faces

	double      				Rho;				    	// density
    double      				R;							// radius of sphere
    double      				Rs;							// radius of sphere factor
	double      				M;					        // mass
	double 						Vol;						// volume
	double      				Kn;					        // normal stiffness
	double      				Kt;					        // tangential stiffness
	double						Gn;							// normal viscous damping coefficient
	double						Gt;							// tangential viscous damping coefficient
	double 						Young;						// Young's modus
	double						Poisson;					// Possion ratio

	Vector3i					Max;						// Max corner of the surrounding box 
	Vector3i 					Min;						// Min corner of the surrounding box 
	Vector3d					BoxL;						// Denmention of the surrounding box
	VectorXd					Coef0;						// Coefficient for 2D Polynomial Particle
	VectorXd					Coef;						// Coefficient for 2D Polynomial Particle

	vector<Vector3d>			P0;				        	// list of point positions at init
	vector<Vector3d>			Ps;				        	// list of point positions at init under spherical coordinate
	vector<Vector3d>			P;				        	// list of point positions at current time step
	vector<Vector2i>			Edges;				        // list of edges
	vector<VectorXi>			Faces;				        // list of faces
	
	vector<Vector3d> 			MetaP0;						// Metaball control points
	vector<Vector3d> 			MetaP;						// Metaball control points
	VectorXd 					MetaK;						// Metaball coefficients
	double 						MetaCC;						// Cricital meta function value for checking cube meta intersection 

	Vector3d					X0;				            // init position
	Vector3d					X;				            // position
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
	Vector3d					Vf;				        	// fixed velocity
	Vector3d					Vc;				        	// constrained velocity
	Vector3d					Wf;				        	// fixed angylar velocity
	Vector3d					Normal;						// normal direction of wall

    bool        				removed;                	// flag for removed DEM_PARTICLEs
	bool				        fixV;				    	// flag for fixed translational velocity
	bool        				fixW;				    	// flag for fixed angular velocity
	bool 						fixed;						// flag for fixed particle with zero velocity
	bool 						crossing[3];
	bool 						crossingFlag;
	bool 						constrained[3];			

    Quaterniond 				Q;				        	// quaternion that describes the rotation
    Quaterniond 				Q0;							// inital quaternion
    Quaterniond					Qf;							// final quaternion (Q0*Q, from object frame to lab frame)
    Quaterniond					Qfi;						// inverse of Qf, used to rotate forces to object frame

    vector< vector<int> >		Lb;							// List of boundary LBM nodes
    vector< vector<int> >		Ln;							// List of neighbours of boundary LBM nodes
    vector< VectorXd >			Lq;							// List of neighbours of boundary LBM nodes

    vector< size_t >			Lp;							// List of particles ID which belong to this group
    vector< double >			Ld;							// List of distance between boundary nodes and particle surFaces for NEBB
    vector< Vector3d >			Li;							// List of position of interpation points for boundary nodes
    vector< Vector3d >			Xmir;				        // mirror positions

    vector< int >				Lcid;
};

inline DEM_PARTICLE::DEM_PARTICLE(int tag, const Vector3d& x, double rho)
{
	PID 	= 0;
	CID 	= 0;
	CIDt	= 0;
    Type	= 0;
	ID		= 0;
	Tag		= tag;
	Group 	= -1;
	MID 	= 0;

	X 		= x;
	X0 		= X;
	Xmir.resize(0);
	Rho		= rho;
	R 		= 0.;
	Rs 		= 0.;
	M 		= 0.;
	Kn	    = 1.0e3;
	Kt		= 2.0e2;
	Gn      = 0.05;
	Gt      = 0.;
	Young 	= 0.;
	Poisson = 0.35;

	Max.setZero();
	Min.setZero();
	BoxL.setZero();

	V.setZero();
	W.setZero();
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
	Vc.setZero();
	Normal.setZero();

	Q.w() = 1;
	Q.vec() << 0.,0.,0.;

	Q0.w() = 1;
	Q0.vec() << 0.,0.,0.;

	fixV	= false;
	fixW	= false;
	removed = false;
	fixed	= false;
	crossing[0] = false;
	crossing[1] = false;
	crossing[2] = false;
	crossingFlag= false;

	constrained[0] = false;
	constrained[1] = false;
	constrained[2] = false;

	Lb.resize(0);
	Ln.resize(0);
	Lq.resize(0);
	Lp.resize(0);
	Lcid.resize(0);

	P0.resize(0);
	Ps.resize(0);
	P.resize(0);
	Edges.resize(0);
	Faces.resize(0);
}

inline void DEM_PARTICLE::FixV(Vector3d& v)
{
	fixV	= true;
	Vf	    = v;
}

inline void DEM_PARTICLE::FixW(Vector3d& w)
{
	fixW	= true;
	Wf	= Q.inverse()._transformVector(w);
}

inline void DEM_PARTICLE::Fix()
{
	fixed = true;
	Vector3d a = Vector3d::Zero();
	FixV(a);
	FixW(a);
}

inline void DEM_PARTICLE::UnFix()
{
	fixed = false;
	fixV = false;
	fixW = false;
}

inline void DEM_PARTICLE::SetG(Vector3d& g)
{
    G       = g;
}

// Velocity Verlet intergrator
// Based on "MBN Explorer Users’ Guide Version 3.0"
// inline void DEM_PARTICLE::VelocityVerlet(double dt, bool movePoints)
// {
// 	//store the position and velocity which before updated
// 	Vector3d Vb, Wb;				//subscript 'b' means before move.
// 	Quaterniond Qb;

// 	Xb	= X;
// 	Qb	= Q;

// 	if (fixV)
// 	{
// 		Vb	= Vf;
// 		Fh.setZero();
// 		Fc.setZero();
// 		G. setZero();
// 		Avb.setZero();
// 	}
// 	else	Vb	= V;	

// 	if (fixW)
// 	{
// 		Wb	= Wf;
// 		Awb.setZero();
// 		Th.setZero();
// 		Tc.setZero();
// 		Awb.setZero();
// 	}
// 	else	Wb	= W;
	
// 	//Update the position and velocity
// 	Vector3d Av	= (Fh + Fc + Fex)/M + G;
// 	// 5.39
// 	X	= Xb + dt*Vb + 0.5*dt*dt*Avb;
// 	// 5.38 and 5.50 
// 	V	= Vb + 0.5*dt*(Avb + Av);

// 	// if (Fc(2)<0.)
// 	// if (Av(2)<-0.01)
// 	if (Fc.norm()>0.)
// 	{
// 		cout << "V: " << V.transpose() << endl;
// 		cout << "Vb: " << Vb.transpose() << endl;
// 		cout << "Avb: " << Avb.transpose() << endl;
// 		cout << "Av: " << Av.transpose() << endl;

// 		cout << "Fh: " << Fh.transpose() << endl;
// 		cout << "Fc: " << Fc.transpose() << endl;
// 		cout << "Fex: " << Fex.transpose() << endl;
// 		cout << "dt: " << dt << endl;
// 		if (Fc(2)<0.)	abort();
// 	}

// 	//Update quaternion
// 	Quaterniond Aq, Qwb, Qawb, dQ, rQ;

// 	Qwb.w()		= 0.;
// 	Qwb.vec()	= 0.5*Wb;
// 	Qawb.w()	= 0.;
// 	Qawb.vec()	= 0.5*Awb;
// 	// 5.43
// 	Aq	= Qb*Qwb;
// 	// 5.44 and 5.46
// 	Q.coeffs()	= Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Aq*Qwb).coeffs() + (Qb*Qawb).coeffs());
// 	// 5.47
// 	Q.normalize();
// 	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
// 	Qfi = Qf.inverse();
// 	// for meta control points
// 	for (size_t i=0; i<MetaP.size(); ++i)
// 	{
// 		MetaP[i] = Qf._transformVector(MetaP0[i]);
// 		MetaP[i] += X;
// 	}
// 	if (movePoints)
// 	{
// 		for (size_t i=0; i<P.size(); ++i)
// 		{
// 			P[i] = Qf._transformVector(P0[i]);
// 			P[i] += Xb;
// 		}
// 	}
// 	//Update the angular velocity
// 	Vector3d Aw0 = I.asDiagonal().inverse()*((Th + Tc + Tex));
// 	// 5.45 and 5.54
// 	Vector3d w0	= Wb + 0.5*dt*(Awb + Aw0);
// 	//First order correction for angular velocity
// 	// 5.55-57
// 	Vector3d Aw1 = I.asDiagonal().inverse()*(-w0.cross(I.asDiagonal()*w0));
// 	// 5.58
// 	Vector3d w1	= w0 + 0.5*dt*Aw1;
// 	//Second order correction for angular velocity
// 	// 5.59-61
// 	Vector3d Aw2 = I.asDiagonal().inverse()*(-w1.cross(I.asDiagonal()*w1));
// 	// 5.62
// 	W	= w1 + 0.5*dt*Aw2;
// 	//store the acceleration for next update
// 	Avb	= Av;
// 	Awb	= Aw0+Aw1+Aw2;

// 	if (constrained[0])		Constrain(0,dt);
// 	if (constrained[1])		Constrain(1,dt);
// 	if (constrained[2])		Constrain(2,dt);
// 	// UpdateBox();
// }

inline void DEM_PARTICLE::UpdatePosition(double dt, bool movePoints, bool resetHydro)
{
	if (fixV)
	{
		V	= Vf;
		Avb.setZero();
	}
	if (fixW)
	{
		W	= Wf;
		Awb.setZero();
	}

	X += dt*V + 0.5*dt*dt*Avb;
	//Update quaternion
	Quaterniond Aq, Qwb, Qawb;

	Qwb.w()		= 0.;
	Qwb.vec()	= 0.5*W;
	Qawb.w()	= 0.;
	Qawb.vec()	= 0.5*Awb;
	Quaterniond Qb = Q;
	// 5.43
	Aq	= Qb*Qwb;
	// 5.44 and 5.46
	Q.coeffs()	= Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Aq*Qwb).coeffs() + (Qb*Qawb).coeffs());
	Q.normalize();
	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();
	// for meta control points
	if (Type==4)
	{
		for (size_t i=0; i<MetaP.size(); ++i)
		{
			MetaP[i] = Qf._transformVector(MetaP0[i]);
			MetaP[i] += X;
		}
	}

	if (movePoints)
	{
		for (size_t i=0; i<P.size(); ++i)
		{
			P[i] = Qf._transformVector(P0[i]);
			P[i] += X;
		}
	}
	// reset contact force
	Fc.setZero();
	Tc.setZero();
	// reset hydro force
	if (resetHydro)
	{
		Fh.setZero();
		Th.setZero();
	}
}

inline void DEM_PARTICLE::UpdateVelocity(double dt)
{
	if (fixV)
	{
		Fh.setZero();
		Fc.setZero();
		Fex.setZero();
		G. setZero();
	}	
	if (fixW)
	{
		Th.setZero();
		Tc.setZero();
		Tex.setZero();
	}
	// Acceleration at current time step
	Vector3d Av	= (Fh + Fc + Fex)/M + G;
	// Update velocity
	V += 0.5*dt*(Avb+Av);
	// if (Tag==4 && Fc.norm()>0.)
	// {
	// 	cout << "M: " << M << endl;
	// 	cout << "fc: " << Fc.transpose() << endl;
	// 	cout << "Avb: " << Avb.transpose() << endl;
	// 	cout << "Av: " << Av.transpose() << endl;
	// 	cout << "V: " << V.transpose() << endl;
	// 	cout << "+++++++++++++++" << endl;
	// }
	// Total torque 
	Vector3d Tt = Th + Tc + Tex;
	//Update the angular velocity 5.51~53
	Vector3d Aw0 = I.asDiagonal().inverse()*Tt;
	// 5.45 and 54
	Vector3d w0	= W + 0.5*dt*(Awb + Aw0);
	//First order correction for angular velocity
	// 5.55-57
	Vector3d Aw1 = I.asDiagonal().inverse()*(-w0.cross(I.asDiagonal()*w0));
	// 5.58
	Vector3d w1	= w0 + 0.5*dt*Aw1;
	//Second order correction for angular velocity
	// 5.59-61
	Vector3d Aw2 = I.asDiagonal().inverse()*(-w1.cross(I.asDiagonal()*w1));
	// 5.62
	W	= w1 + 0.5*dt*Aw2;
	// Store information for next update
	Avb	= Av;
	Awb = I.asDiagonal().inverse()*Tt;
}

inline void DEM_PARTICLE::Constrain(int d, double dt)
{
	X(d) = Xb(d)+Vc(d)*dt;
	V(d) = Vc(d);
}

inline void DEM_PARTICLE::UpdateBox(size_t D)
{
	//Update the range of wrapping box
	for (size_t d=0; d<D; ++d)
	{
		Max(d) = (int) (X(d)+BoxL(d));
		Min(d) = (int) (X(d)-BoxL(d));
	}
	// Xmir.clear();
	crossing[0] = false;
	crossing[1] = false;
	crossing[2] = false;
	crossingFlag= false;
}

inline void DEM_PARTICLE::ZeroForceTorque(bool h, bool c)
{
	//set zero hydro and contact force and torque after moved.
	if (h)
	{
		Fh.setZero();
		Th.setZero();
	}
	if (c)
	{
		Fc.setZero();
		Tc.setZero();
	}
}

inline double DEM_PARTICLE::CalEnergy()
{
	double eng = 0.5*(M*V.squaredNorm()+I(0)*W(0)*W(0)+I(1)*W(1)*W(1)+I(2)*W(2)*W(2));
	return eng;
}

inline Vector3d DEM_PARTICLE::CalAngularMoment()
{
	Vector3d am = M*X.cross(V);
	return am;
}

inline void DEM_PARTICLE::RerangeFaceElementOrder() // only works for convex shape
{
	// rerange face elements order (right hand role for outside normal)
	for (size_t i=0; i<Faces.size();++i)
	{
		size_t ind;	// index of the point that not belongs to this face i
		for (size_t j=0; j<P.size(); ++j)
		{
			bool exist = false;
			for (int k=0; k<Faces[i].size();++k)
			{
				if (j-Faces[i](k)==0)	// check if j belongs to face i
				{
					exist = true;
					break;
				}
			}
			if (!exist)
			{
				ind = j;
				break;
			}
		}
		Vector3d v0 = P[Faces[i](0)]-P[Faces[i](1)];
		Vector3d v1 = P[Faces[i](2)]-P[Faces[i](1)];
		Vector3d v3 = P[ind]-P[Faces[i](1)];
		if ((v1.cross(v0)).dot(v3)>0)	// check norm is with same side of point ind
		{
			Vector3i facer = Faces[i].reverse();	// if so reverse order of face i
			Faces[i] = facer;
		}
	}
}

inline void DEM_PARTICLE::SetPlane(Vector3d& n)
{
	Type = -1;
	Normal = n;
	M = 1.0e12;
	I(0) = I(1) = I(2) = 1.0e12;
}

inline void DEM_PARTICLE::SetSphere(double r)
{
    Type= 1;
	R	= r;
	Vol = 4./3.*M_PI*R*R*R;
	M	= Rho*Vol;
	I(0)	= 0.4*M*R*R;
	I(1)	= I(0);
	I(2)	= I(0);

	P0.resize(1);
	P0[0] << R, 0., 0.;
	P.resize(1);
	P[0] << R, 0., 0.;

	BoxL << R+1, R+1, R+1;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));

	MetaP.push_back(X);
	MetaK.resize(1);
	MetaK(0) = R*R;
}
// only used for 2d
inline void DEM_PARTICLE::SetDisk2D(double r)
{
    Type= 2;
	R	= r;
	Vol = M_PI*R*R;
	M	= Rho*Vol;
	I(0)	= 0.25*M*R*R;
	I(1)	= I(0);
	I(2)	= 2.*I(0);

	P0.resize(1);
	P0[0] << R, 0., 0.;
	P.resize(1);
	P[0] << R, 0., 0.;

	BoxL << R+1, R+1, 0;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
}

inline void DEM_PARTICLE::SetCuboid(double lx, double ly, double lz)
{
    Type= 3;
	R	= 0.1;
	Vol = lx*ly*lz;
	M	= Rho*Vol;
	I(0)	= M*(ly*ly+lz*lz)/12.;
	I(1)	= M*(lx*lx+lz*lz)/12.;
	I(2)	= M*(lx*lx+ly*ly)/12.;

	P0.resize(8);
	P0[0] << -0.5*lx, -0.5*ly, -0.5*lz;
	P0[1] <<  0.5*lx, -0.5*ly, -0.5*lz;
	P0[2] <<  0.5*lx,  0.5*ly, -0.5*lz;
	P0[3] << -0.5*lx,  0.5*ly, -0.5*lz;
	P0[4] << -0.5*lx, -0.5*ly,  0.5*lz;
	P0[5] <<  0.5*lx, -0.5*ly,  0.5*lz;
	P0[6] <<  0.5*lx,  0.5*ly,  0.5*lz;
	P0[7] << -0.5*lx,  0.5*ly,  0.5*lz;
	P.resize(P0.size());
	for (size_t i=0; i<P.size(); ++i)
	{
		P[i] = P0[i]+X;
	}

	VectorXi face(4);
	face << 0, 3, 2, 1;
	Faces.push_back(face);
	face << 4, 5, 6, 7;
	Faces.push_back(face);
	face << 1, 5, 4, 0;
	Faces.push_back(face);
	face << 2, 3, 7, 6;
	Faces.push_back(face);
	face << 1, 2, 6, 5;
	Faces.push_back(face);
	face << 0, 4, 7, 3;
	Faces.push_back(face);

	Nfe = 30;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetTetrahedron(vector<Vector3d> ver)
{
    Type= 3;
	R	= 0.1;
	Vol	= abs((ver[0]-ver[3]).dot((ver[1]-ver[3]).cross(ver[2]-ver[3])))/6.;
	M	= Rho*Vol;
	X 	= (ver[0]+ver[1]+ver[2]+ver[3])/4.;
	P = ver;
	// add Faces
	VectorXi face(3);
	face << 0, 1, 2;
	Faces.push_back(face);
	face << 1, 2, 3;
	Faces.push_back(face);
	face << 2, 3, 0;
	Faces.push_back(face);
	face << 0, 3, 1;
	Faces.push_back(face);
	RerangeFaceElementOrder();

	Matrix3d Ip;	// Inertia tensor
	CalVolumeIntegrals(Faces, P, Rho, Vol, M, X, Ip);
	// move P0 to make sure mass centre is the origin
	P0.resize(P.size());
	for (size_t i=0; i<P.size(); ++i)	P0[i] = P[i]-X;
	Vector3d xc;
	CalVolumeIntegrals(Faces, P0, Rho, Vol, M, xc, Ip);
	// use eigen vectors to find the rotation from lab frame to object frame
	SelfAdjointEigenSolver<Matrix3d> sol (Ip);
	Matrix3d em = sol.eigenvectors();
	if ((em.col(0).cross(em.col(1))).dot(em.col(2))<0)	em.col(2) *= -1;
	Quaterniond q0(em);
	// rotation P0 to object frame
	for (size_t i=0; i<P.size(); ++i)	P0[i] = q0.inverse()._transformVector(P0[i]);

	Q0 = q0;
	Nfe = 16;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetMetaball(double rs, vector<Vector3d>& metaP, VectorXd& metaK, vector<Vector3d>& unitSphereP, vector<VectorXi>& unitSphereF)
{
	Type = 4;
	Rs = rs;
	MetaP = metaP;
	MetaK = metaK;
	Faces = unitSphereF;
	GetMetaSurfaceMesh(unitSphereP, MetaP, MetaK, P);
	Matrix3d Ip;	// Inertia tensor
	CalVolumeIntegrals(Faces, P, Rho, Vol, M, X, Ip);
	// use eigen vectors to find the rotation from lab frame to object frame
	SelfAdjointEigenSolver<Matrix3d> sol (Ip);
	Matrix3d em = sol.eigenvectors();
	I = sol.eigenvalues();
	if ((em.col(0).cross(em.col(1))).dot(em.col(2))<0)	em.col(2) *= -1;
	Quaterniond q0(em);
	// move P0 to make sure mass centre is the origin and rotation P0 to object frame
	P0.resize(P.size());
	double maxr = 0.;
	double minc = 1.e12;
	for (size_t i=0; i<P.size(); ++i)
	{
		P0[i] = P[i]-X;
		P0[i] = q0.inverse()._transformVector(P0[i]);
		double r = P0[i].norm();
		if (r>maxr)	maxr = r;
		Vector3d pm = P0[i]+(0.87+Rs)*P0[i].normalized();
		double c = CalMetaC(metaP, metaK, pm);
		if (c<minc)	minc = c;
	}
	R = 1.1*(maxr+Rs);
	MetaCC = minc;
	MetaP0.resize(MetaP.size());
	for (size_t i=0; i<MetaP.size(); ++i)
	{
		MetaP0[i] = MetaP[i]-X;
		MetaP0[i] = q0.inverse()._transformVector(MetaP0[i]);
	}
	Q0 = q0;
	Nfe = 4*unitSphereF.size();

	BoxL << R, R, R;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
	// cout << "MetaCC: " << MetaCC << endl;
	// abort();
}

inline void DEM_PARTICLE::SetDrum(double rd, Vector3d& xd)
{
	Type= -2;
	R	= rd;
	X 	= xd;
}