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
	void SetSphere(double r);								// Change DEM_PARTICLE to Sphere
	void SetDisk2D(double r);								// Change DEM_PARTICLE to 2D disk
	void Set2DPolynomialParticle(VectorXd& coef);			// Change DEM_PARTICLE to 2D Polynomial Particle
	void SetCuboid(double lx, double ly, double lz);		// Change DEM_PARTICLE to Cuboid
	void UpdateCoef();
	// void DistanceToSurface(Vector3d& x);

    int         				Type;                       // Type of DEM_PARTICLE, for 0 is disk2d, for 1 is sphere or 2 is cube etc.
	int         				ID; 				    	// index of DEM_PARTICLE in the list 
	int         				Tag;				    	// tag of DEM_PARTICLE
	int         				Group;				    	// tag of Group
	int 						MID;						// material ID is used to find friction coefficient 

	double      				Rho;				    	// density
    double      				R;							// radius of sphere
	double      				M;					        // mass
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
	Vector3d					Th;				        	// Hydro torque under lab frame
	Vector3d					Tc;				        	// Contact torque under lab frame
	Vector3d					Tex;			        	// Variable external torque under lab frame
	Vector3d					Avb;				        // acceleration of velocity before
	Vector3d					Awb;				        // acceleration of angluar velocity before under DEM_PARTICLE frame
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
	Quaterniond Aq, Qwb, Qawb, deltaQ, rQ;

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

	//Update vertices
	// deltaQ	= Qb.inverse()*Q;

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

inline void DEM_PARTICLE::SetSphere(double r)
{
    Type= 1;
	R	= r;
	M	= 4./3.*Rho*M_PI*R*R*R;
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
    Type= 1;
	R	= r;
	M	= Rho*M_PI*R*R;
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

// inline double DEM_PARTICLE::DistanceToSurface(Vector3d& x)
// {
// 	double dis = (X-x).norm()-R;
// 	for (size_t i=0; i<Xmir.size(); ++i)
// 	{
// 		dis = min(dis, (Xmir[i]-x).norm()-R);
// 	}
// 	return dis;
// }
// inline void DEM_PARTICLE::Set2DPolynomialParticle(VectorXd& coef)
// {
//     Type= 4;
//     Coef0 = coef;
// 	Coef = coef;
// 	// if (Coef.size()!=15)
// 	// {
// 	// 	cout << "Wrong size of Coefficient for PDEM" << endl;
// 	// 	abort();
// 	// }
// 	R	= 1.7;
// 	M	= 4./3.*Rho*M_PI*R*R*R;
// 	I(0)	= 0.4*M*R*R;
// 	I(1)	= I(0);
// 	I(2)	= I(0);

// 	size_t np = 20;

// 	P0.resize(np+1);
// 	P.resize(np+1);
// 	// VectorXi face(np);
// 	P0[0] << 0.724494, 0., 0.;
// 	P0[1] << 0.634993, 0.206322, 0.;
// 	P0[2] << 0.512132, 0.372086, 0.;
// 	P0[3] << 0.370052, 0.508333, 0.;
// 	P0[4] << 0.207438, 0.63843, 0.;
// 	P0[5] << 0., 0.786151, 0.;
// 	P0[6] << -0.321996, 0.991001, 0.;
// 	P0[7] << -0.882223, 1.21428, 0.;
// 	P0[8] << -1.31978, 0.958876, 0.;
// 	P0[9] << -1.30216, 0.423098, 0.;
// 	P0[10] << -1.22074, 0., 0.;
// 	P0[11] << -1.1452, -0.372099, 0.;
// 	P0[12] << -0.99897, -0.725794, 0.;
// 	P0[13] << -0.64472, -0.887381, 0.;
// 	P0[14] << -0.275359, -0.847468, 0.;
// 	P0[15] << -0., -0.786151, 0.;
// 	P0[16] << 0.236186, -0.726905, 0.;
// 	P0[17] << 0.468913, -0.645403, 0.;
// 	P0[18] << 0.670926, -0.487456, 0.;
// 	P0[19] << 0.753471, -0.244818, 0.;
// 	P0[20] << 0., -0., 0.;

// 	for (size_t i=0; i<np; ++i)
// 	{
// 		P[i] = P0[i]+X;
// 		Vector3i face (20, i, (i+1)%np);
// 		Faces.push_back(face);
// 		// if (i==1)	Faces.push_back(face);
// 	}
// 	P[20] = P0[20]+X;

// 	BoxL << R+1, R+1, R+1;

// 	Max(0) 	= (int) (X(0)+BoxL(0));
// 	Max(1) 	= (int) (X(1)+BoxL(1));

// 	Min(0) 	= (int) (X(0)-BoxL(0));
// 	Min(1) 	= (int) (X(1)-BoxL(1));
// }

// inline void DEM_PARTICLE::UpdateCoef()
// {
// 	// Displace
// 	double u = X(0)-X0(0);
// 	double v = X(1)-X0(1);
// 	double u2 = u*u;
// 	double uv = u*v;
// 	double v2 = v*v;

// 	VectorXd t = coef;

// 	coef(5) = t(5)-4.*t(0)*u-t(1)*v;
// 	coef(6) = -3.*t(1)*u-2.*t(2)*v+t(6);
// 	coef(7) = -2.*t(2)*u-3.*t(3)*v+t(7);
// 	coef(8) = -t(3)*u-4.*t(4)*v+t(8);
// 	coef(9) = 6.*t(0)*u2+3.*t(1)*uv-3.*t(5)*u+t(2)*v2-t(6)*v+t(9);
// 	coef(10) = 3.*t(1)*u2+4.*t(2)*uv-2.*t(6)*u+3.*t(3)*v2-2.*t(7)*v+t(10);
// 	coef(11) = t(2)*u2+3.*t(3)*uv-t(7)*u+6.*t(4)*v2-3.*t(8)*v+t(11);
// 	coef(12) = -4.*t(0)*u*u2-3.*t(1)*u2*v+3.*t(5)*u2-2.*t(2)*u*v2+2.*t(6)*uv-2.*t(9)*u-t(3)*v2*v+t(7)*v2-t(10)*v+t(12);
// 	coef(13) = -t(1)*u2*u-2*t(2)*u2*v+t(6)*u2-3*t(3)*u*v2+2*t(7)*uv-t(10)*u-4*t(4)*v2*v+3*t(8)*v2-2*t(11)*v+t(13);
// 	coef(14) =  t(0)*u2*u2+t(1)*u2*uv-t(5)*u2*u+t(2)*u2*v2-t(6)*u2*v+t(9)*u2+t(3)*uv*v2-t(7)*u*v2+t(10)*uv-t(12)*u+t(4)*v2*v2-t(8)*v2*v+t(11)*v2-t(13)*v+t(14);
// }

// inline void DEM_PARTICLE::UpdateCoef()
// {
// 	// Displace
// 	double u = X(0)-X0(0);
// 	double v = X(1)-X0(1);
// 	Coef.setZero();
// 	Coef(0) = 1.;
// 	Coef(4) = 1.;
// 	Coef(5) = -4.*u;
// 	Coef(8) = -4.*v;
// 	Coef(9) = 6.*u*u;
// 	Coef(11) = 6.*v*v-u+1.;
// 	Coef(12) = v*v-v-4.*u*u*u+1.;
// 	Coef(13) = -4.*v*v*v+2.*u*v-u-2.*v;
// 	Coef(7) = 1.;
// 	Coef(10) = -2.*v+1.;
// 	Coef(14) = u*u*u*u+v*v*v*v-v*v*u+u*v+v*v-u-1.;
// }
// inline void DEM_PARTICLE::SetCuboid(double lx, double ly, double lz)
// {
//     Type= 2;
// 	R	= 0.5*sqrt(lx*lx+ly*ly+lz*lz);
// 	M	= Rho*lx*ly*lz;
// 	I(0)	= M*(ly*ly+lz*lz);
// 	I(1)	= M*(lx*lx+lz*lz);
// 	I(2)	= M*(lx*lx+ly*ly);

// 	P0.resize(8);
// 	P0[0] << -0.5*lx, -0.5*ly, -0.5*lz;
// 	P0[1] <<  0.5*lx, -0.5*ly, -0.5*lz;
// 	P0[2] <<  0.5*lx,  0.5*ly, -0.5*lz;
// 	P0[3] << -0.5*lx,  0.5*ly, -0.5*lz;
// 	P0[4] << -0.5*lx, -0.5*ly,  0.5*lz;
// 	P0[5] <<  0.5*lx, -0.5*ly,  0.5*lz;
// 	P0[6] <<  0.5*lx,  0.5*ly,  0.5*lz;
// 	P0[7] << -0.5*lx,  0.5*ly,  0.5*lz;
// 	P.resize(P0.size());
// 	for (size_t i=0; i<P.size(); ++i)
// 	{
// 		P[i] = P0[i]+X;
// 	}

// 	VectorXi face(4);
// 	face << 0, 1, 2, 3;
// 	Faces.push_back(face);
// 	face << 4, 5, 6, 7;
// 	Faces.push_back(face);
// 	face << 0, 1, 5, 4;
// 	Faces.push_back(face);
// 	face << 2, 3, 7, 6;
// 	Faces.push_back(face);
// 	face << 1, 2, 6, 5;
// 	Faces.push_back(face);
// 	face << 0, 3, 7, 4;
// 	Faces.push_back(face);

// 	Max(0) 	= (int) (X(0)+BoxL(0));
// 	Max(1) 	= (int) (X(1)+BoxL(1));
// 	Max(2) 	= (int) (X(2)+BoxL(2));

// 	Min(0) 	= (int) (X(0)-BoxL(0));
// 	Min(1) 	= (int) (X(1)-BoxL(1));
// 	Min(2) 	= (int) (X(2)-BoxL(2));
// }

// inline void DEM_PARTICLE::SetPlane(double r)
// {
//     Type= 1;
// 	R	= r;
// 	M	= 4./3.*Rho*M_PI*R*R*R;
// 	I(0)	= 0.4*M*R*R;
// 	I(1)	= I(0);
// 	I(2)	= I(0);

// 	P0.resize(1);
// 	P0[0] << R, 0., 0.;
// 	P.resize(1);
// 	P[0] << R, 0., 0.;

// 	MaxX	= (int) (X(0)+R+1);
// 	MaxY	= (int) (X(1)+R+1);
// 	MaxZ	= (int) (X(2)+R+1);

// 	MinX	= (int) (X(0)-R);
// 	MinY	= (int) (X(1)-R);
// 	MinZ	= (int) (X(2)-R);
// }

/*inline void DEM_PARTICLE::SetPolyhedron(double lx, double ly, double lz)
{
    Type= 3;
	R	= 0.5*sqrt(lx*lx+ly*ly+lz*lz);
	M	= Rho*lx*ly*lz;
	I(0)	= M*(ly*ly+lz*lz);
	I(1)	= M*(lx*lx+lz*lz);
	I(2)	= M*(lx*lx+ly*ly);

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
		P[i] = P0[i];
	}

	MaxX	= (int) (X(0)+R+1);
	MaxY	= (int) (X(1)+R+1);
	MaxZ	= (int) (X(2)+R+1);

	MinX	= (int) (X(0)-R);
	MinY	= (int) (X(1)-R);
	MinZ	= (int) (X(2)-R);
}*/

// class DISK2D:public DEM_PARTICLE
// {
// public:
// 	DISK2D(int tag, double r, const Vector3d& x, double rho);
// };

// inline DISK2D::DISK2D(int tag, double r, const Vector3d& x, double rho):DEM_PARTICLE(tag, x, rho)
// {
//     Type= 0;
// 	R	= r;
// 	M	= Rho*M_PI*R*R;
// 	I(0)	= 0.25*M*R*R;
// 	I(1)	= I(0);
// 	I(2)	= 2.*I(0);

// 	P0.resize(1);
// 	P0[0] << 1., 0., 0.;
// 	P.resize(1);
// 	P[0] << 1., 0., 0.;

// 	MaxX	= (int) (X(0)+R+1);
// 	MaxY	= (int) (X(1)+R+1);
// 	MaxZ	= 0;

// 	MinX	= (int) (X(0)-R);
// 	MinY	= (int) (X(1)-R);
// 	MinZ	= 0;
// }