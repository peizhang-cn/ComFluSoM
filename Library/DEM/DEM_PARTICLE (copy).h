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

class DEM_PARTICLE						    				// class for a single DEM_PARTICLE
{
public:
	DEM_PARTICLE(int tag, const Vector3d& x, double rho);
	void Set(double kn);			    					// set physical parameters of the DEM_PARTICLEs
	void SetG(Vector3d& g);			    					// set external acceleration
	void VelocityVerlet(double dt);							// move the DEM_PARTICLE based on Velocity Verlet intergrator
	void UpdateBox(size_t D);								// Update the range of warpping box and reset cross flags
	void Constrain(int d, double dt);						// Constrain on d direction with v velocity, applied after VelocityVerlet
	void FixV(Vector3d& v);			    					// fix the velocity
	void FixW(Vector3d& w);			    					// fix the angular velocity
	void Fix();												// fix particle with zero velocity
	void UnFix();											// reset all fixs flags
	void ZeroForceTorque(bool h, bool c);			    	// set zero force and torque
	void RerangeFaceElementOrder();							// rerange face element order to make sure norm towards to outside, only works for convex shapes
	void SetSphere(double r);								// Change DEM_PARTICLE to Sphere
	void SetDisk2D(double r);								// Change DEM_PARTICLE to 2D disk
	void Set2DPolynomialParticle(VectorXd& coef);			// Change DEM_PARTICLE to 2D Polynomial Particle
	void SetCuboid(double lx, double ly, double lz);		// Change DEM_PARTICLE to Cuboid
	void SetTetrahedron(vector<Vector3d> ver);
	void UpdateCoef();
	// void DistanceToSurface(Vector3d& x);

    int         				Type;                       // Type of DEM_PARTICLE, for 0 is disk2d, for 1 is sphere or 2 is cube etc.
	int         				ID; 				    	// index of DEM_PARTICLE in the list 
	int         				Tag;				    	// tag of DEM_PARTICLE
	int         				Group;				    	// tag of Group
	int 						MID;						// material ID is used to find friction coefficient
	int 						Nfe;						// Total number of elements in faces

	double      				Rho;				    	// density
    double      				R;							// radius of sphere
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

    Quaterniond 				Q;				        	// quaternion that describes the orientation of DEM_PARTICLE frame inrespect to the lab frame
    Quaterniond 				Q0;

    vector< vector<int> >		Lb;							// List of boundary LBM nodes
    vector< vector<int> >		Ln;							// List of neighbours of boundary LBM nodes
    vector< VectorXd >			Lq;							// List of neighbours of boundary LBM nodes

    vector< size_t >			Lp;							// List of particles ID which belong to this group
    vector< double >			Ld;							// List of distance between boundary nodes and particle surFaces for NEBB
    vector< Vector3d >			Li;							// List of position of interpation points for boundary nodes
    vector< Vector3d >			Xmir;				        // mirror positions
};

inline DEM_PARTICLE::DEM_PARTICLE(int tag, const Vector3d& x, double rho)
{
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
	M 		= 0.;
	Kn	    = 1.0e3;
	Kt		= 2.0e2;
	Gn      = 0.05;
	Gt      = 0.;
	Young 	= 0.;
	Poisson = 0.;

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
inline void DEM_PARTICLE::VelocityVerlet(double dt)
{
	//store the position and velocity which before updated
	Vector3d Vb, Wb;				//subscript 'b' means before move.
	Quaterniond Qb;

	Xb	= X;
	Qb	= Q;

	if (fixV)
	{
		Vb	= Vf;
		Fh.setZero();
		Fc.setZero();
		G. setZero();
		Avb.setZero();
	}
	else	Vb	= V;	

	if (fixW)
	{
		Wb	= Wf;
		Awb.setZero();
		Th.setZero();
		Tc.setZero();
		Awb.setZero();
	}
	else	Wb	= W;
	
	//Update the position and velocity
	Vector3d Av	= (Fh + Fc + Fex)/M + G;
	// 5.39
	X	= Xb + dt*Vb + 0.5*dt*dt*Avb;
	// 5.38 and 5.50 
	V	= Vb + 0.5*dt*(Avb + Av);

	//Update quaternion
	Quaterniond Aq, Qwb, Qawb, dQ, rQ;

	Qwb.w()		= 0.;
	Qwb.vec()	= 0.5*Wb;
	Qawb.w()	= 0.;
	Qawb.vec()	= 0.5*Awb;
	// 5.43
	Aq	= Qb*Qwb;
	// 5.44 and 5.46
	Q.coeffs()	= Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Aq*Qwb).coeffs() + (Qb*Qawb).coeffs());
	// 5.47
	Q.normalize();

	// dQ.w() = cos(Wf(0)/2.);
	// dQ.vec() << 1,0,0;
	// dQ.vec() *= sin(Wf(0)/2.);
	// dQ.normalize();
	// Q = dQ*Q;
	// cout << Q.vec().normalized().transpose() << endl;

	//Update vertices
	// dQ	= Qb.inverse()*Q;

	// for (size_t i=0; i<P.size(); ++i)
	// {
	// 	P[i] = dQ._transformVector(P[i]);
	// }

	for (size_t i=0; i<P.size(); ++i)
	{
		P[i] = Q._transformVector(P0[i]);
		P[i] += X;
	}
	//Update the angular velocity
	Vector3d Aw0 = I.asDiagonal().inverse()*((Th + Tc + Tex));
	// 5.45 and 5.54
	Vector3d w0	= Wb + 0.5*dt*(Awb + Aw0);
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

	//store the acceleration for next update
	Avb	= Av;
	Awb	= Aw0+Aw1+Aw2;

	if (constrained[0])		Constrain(0,dt);
	if (constrained[1])		Constrain(1,dt);
	if (constrained[2])		Constrain(2,dt);
	// UpdateBox();
}

inline void DEM_PARTICLE::Constrain(int d, double dt)
{
	X(d) = Xb(d)+Vc(d)*dt;
	V(d) = Vc(d);
}

inline void DEM_PARTICLE::UpdateBox(size_t D)
{
	//Update the range of warpping box
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

inline void DEM_PARTICLE::RerangeFaceElementOrder() // only works for convex shape
{
	// rerange face elements order (right hand role for outside normal)
	for (size_t i=0; i<Faces.size();++i)
	{
		size_t ind;	// index of the point that not belongs to this face i
		for (size_t j=0; j<P.size(); ++j)
		{
			bool exist = false;
			for (size_t k=0; k<Faces[i].size();++k)
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
	R	= 0.5*sqrt(lx*lx+ly*ly+lz*lz);
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
    Type= 4;
	// R	= 0.5*sqrt(lx*lx+ly*ly+lz*lz);
	Vol	= abs((ver[0]-ver[3]).dot((ver[1]-ver[3]).cross(ver[2]-ver[3])))/6.;
	M	= Rho*Vol;
	X 	= (ver[0]+ver[1]+ver[2]+ver[3])/4.;
	P = ver;

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
	cout << "xc: " << X.transpose() << endl;
	// move P0 to make sure mass centre is the origin
	P0.resize(P.size());
	for (size_t i=0; i<P.size(); ++i)	P0[i] = P[i]-X;
	// use eigen vectors to find the rotation from lab frame to object frame
	SelfAdjointEigenSolver<Matrix3d> sol (Ip);
	Matrix3d em = sol.eigenvectors();
	if ((em.col(0).cross(em.col(1))).dot(em.col(2))<0)	em.col(2) *= -1;
	Quaterniond q0(em);

	cout << em << endl;
	cout << sol.eigenvalues().transpose() << endl;
	cout << em.col(0).dot(em.col(2)) << endl;
	cout << em.col(0).dot(em.col(1)) << endl;
	cout << em.col(1).dot(em.col(2)) << endl;
// 	// abort();

	Vector3d vx0 = em.col(0);
	Vector3d vy0 = em.col(1);
	Vector3d vz0 = em.col(2);

	Vector3d vx (1,0,0);
	Vector3d vy (0,1,0);
	Vector3d vz (0,0,1);
	Quaterniond qx, qy, qz;
	qx=Quaterniond::FromTwoVectors(vx0, vx);
	vy0 = qx._transformVector(vy0);
	// vz0 = qx._transformVector(vz0);
	qy=Quaterniond::FromTwoVectors(vy0, vy);
	vz0 = (qy*qx)._transformVector(vz0);
	cout << "vz0= " << vz0.transpose() <<endl;
	// abort();
	// qz=Quaterniond::FromTwoVectors(em.col(2), vz);

// 	qx=Quaterniond::FromTwoVectors(vx, vy);
// 	cout << "vec: " << qx.vec().normalized().transpose() << endl;
// cout << "-------------" << endl;
// cout << qx.w() << endl;

	cout << "-------------" << endl;
	// cout << qx.coeffs().transpose() << endl;
	// cout << qx.norm() << endl;
	// cout << qy.coeffs().transpose() << endl;
	// cout << qy.norm() << endl;
	// cout << qz.coeffs().transpose() << endl;
	// cout << qz.norm() << endl;

	cout << (qy*qx).coeffs().transpose() << endl;
	// cout << (qy*qx*qz).coeffs().transpose() << endl;
	// cout << (qz*qx*qy).coeffs().transpose() << endl;
	cout << q0.coeffs().transpose() << endl;
	cout << q0.inverse().coeffs().transpose() << endl;

	P0.resize(P.size());
	for (size_t i=0; i<P.size(); ++i)
	{
		// P0[i] = P[i]-X;
		P0[i] = q0.inverse()._transformVector(P0[i]);
		P[i] = P0[i];
		// P[i] -= X;
		// cout << "p0: " << P0[i].transpose() << endl;
	}

	CalVolumeIntegrals(Faces, P, Rho, Vol, M, X, Ip);
	cout << "xc: " << X.transpose() << endl;
	cout <<"122222" << endl;
	// abort();

	// /*Q = */Q0 = q0.inverse();

	// abort();

	Nfe = 16;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}