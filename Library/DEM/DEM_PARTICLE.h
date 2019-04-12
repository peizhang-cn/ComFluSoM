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
	void VelocityVerlet(bool firstStep, double dt);			// move the DEM_PARTICLE based on Velocity Verlet intergrator
	void FixV(Vector3d& v);			    					// fix the velocity
	void FixW(Vector3d& w);			    					// fix the angular velocity
	void ZeroForceTorque();			    					// set zero force and torque

    int         				Type;                       // Type of DEM_PARTICLE, for 0 is disk2d, for 1 is sphere or 2 is cube etc.
	int         				ID; 				    	// index of DEM_PARTICLE in the list 
	int         				Tag;				    	// tag of DEM_PARTICLE

	double      				Rho;				    	// density
    double      				R;							// radius of sphere
	double      				M;					        // mass
	double      				Kn;					        // normal stiffness
	double      				Kt;					        // tangential stiffness
	double						Gn;							// normal viscous coefficient
	double						Gt;							// tangential viscous coefficient
	double						Mu;							// friction coefficient

	int							MaxX;						// Max x position of warpping box for DEM_PARTICLEs
	int							MaxY;						// Max y position of warpping box for DEM_PARTICLEs
	int							MaxZ;						// Max z position of warpping box for DEM_PARTICLEs

	int							MinX;						// Min x position of warpping box for DEM_PARTICLEs
	int							MinY;						// Min y position of warpping box for DEM_PARTICLEs
	int							MinZ;						// Min z position of warpping box for DEM_PARTICLEs

	int							N;							// angle of increasement under spherical coordinates
	ArrayXXd					Ps0;
	vector<Vector3d>			P0;				        	// list of points at init
	vector<Vector3d>			P;				        	// list of points at current time step
	Vector3d					X;				            // position
	Vector3d					Xb;				            // position before move, only used for DELBM to find refilling LBM nodes
	Vector3d					V;				            // velocity in the center
	Vector3d					W;				            // angular velocity under DEM_PARTICLE frame
	Vector3d					I;				            // inertia under DEM_PARTICLE frame
	Vector3d					G;				        	// gravity
	Vector3d					Fh;				        	// hydro force
	Vector3d					Fc;				        	// contact force
	Vector3d					Fex;				        // external force
	Vector3d					Th;				        	// hydro torque under lab frame
	Vector3d					Tc;				        	// contact torque under lab frame
	Vector3d					Tex;			        	// external torque under lab frame
	Vector3d					Avb;				        // acceleration of velocity before
	Vector3d					Awb;				        // acceleration of angluar velocity before under DEM_PARTICLE frame
	Vector3d					Vf;				        	// fixed velocity
	Vector3d					Wf;				        	// fixed angylar velocity

    bool        				removed;                	// flag for removed DEM_PARTICLEs
	bool				        fixV;				    	// flag for fixed translational velocity
	bool        				fixW;				    	// flag for fixed angular velocity

    Quaterniond 				Q;				        	// quaternion that describes the orientation of DEM_PARTICLE frame inrespect to the lab frame
    Quaterniond 				Q0;

    vector< vector<int> >		Lb;							//List of boundary LBM nodes
    vector< vector<int> >		Ln;							//List of neighbours of boundary LBM nodes
    vector< VectorXd >			Lq;							//List of neighbours of boundary LBM nodes

    vector< double >			Ld;							//List of distance between boundary nodes and particle surface for NEBB
    vector< Vector3d >			Li;							//List of position of interpation points for boundary nodes
};

inline DEM_PARTICLE::DEM_PARTICLE(int tag, const Vector3d& x, double rho)
{
    Type	= 0;
	ID		= 0;
	Tag		= tag;

	X 		= x;
	Rho		= rho;
	R 		= 0.;
	M 		= 0.;
	Kn	    = 1.0e3;
	Kt		= 5.;
	Gn      = 0.034;
	Gt      = 0.034;
	Mu 		= 0.15;

	MaxX	= 0.;
	MaxY	= 0.;
	MaxZ	= 0.;
	MinX	= 0.;
	MinY	= 0.;
	MinZ	= 0.;

	V.setZero();
	W.setZero();
	G.setZero();
	I.setZero();
	Fh.setZero();
	Fc.setZero();
	Fex.setZero();
	Th.setZero();
	Tc.setZero();
	Tex.setZero();
	Avb.setZero();
	Awb.setZero();

	Q.w() = 1;
	Q.vec() << 0.,0.,0.;

	Q0.w() = 1;
	Q0.vec() << 0.,0.,0.;

	fixV	= false;
	fixW	= false;
	removed = false;

	Lb.resize(0);
	Ln.resize(0);
	Lq.resize(0);
}

inline void DEM_PARTICLE::FixV(Vector3d& v)
{
	fixV	= true;
	Vf	    = v;
}

inline void DEM_PARTICLE::FixW(Vector3d& w)
{
	fixW	= true;
	Wf	= w;
}

inline void DEM_PARTICLE::SetG(Vector3d& g)
{
    G       = g;
}

// Velocity Verlet intergrator
inline void DEM_PARTICLE::VelocityVerlet(bool firstStep, double dt)
{
	//store the position and velocity which before updated
	Vector3d Xb, Vb, Wb;				//subscript 'b' means before move.
	Quaterniond Qb;

	Xb	= X;
	Qb	= Q;

	if (fixV)
	{
		Vb	= Vf;
		Fh.setZero();
		Fc.setZero();
		G. setZero();
	}
	else	Vb	= V;	

	if (fixW)
	{
		Wb	= Wf;
		Awb.setZero();
		Th.setZero();
		Tc.setZero();
	}
	else	Wb	= W;

	if (firstStep)
	{
		Avb = (Fh + Fc + Fex)/M + G;
		Awb = I.asDiagonal().inverse()*((Th + Tc + Tex));
	}

	//Update the position and velocity
	Vector3d Av	= (Fh + Fc + Fex)/M + G;

	X	= Xb + dt*Vb + 0.5*dt*dt*Avb;				
	V	= Vb + 0.5*dt*(Avb + Av);

	//Update quaternion
	Quaterniond Aq, Qwb, Qawb, deltaQ, rQ;

	Qwb.w()		= 0.;
	Qwb.vec()	= 0.5*Wb;
	Qawb.w()	= 0.;
	Qawb.vec()	= 0.5*Awb;

	Aq	= Qb*Qwb;

	Q.coeffs()	= Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Aq*Qwb).coeffs() + (Qb*Qawb).coeffs());

	Q.normalize();

	//Update vertices
	deltaQ	= Qb.inverse()*Q;

	for (size_t i=0; i<P.size(); ++i)
	{
		P[i] = Q._transformVector(P0[i]);
		P[i] += X;
	}
	//Update the angular velocity
	Vector3d Aw0 = I.asDiagonal().inverse()*((Th + Tc + Tex));
	Vector3d w0	= Wb + 0.5*dt*(Awb + Aw0);
	//First order correction for angular velocity
	Vector3d Aw1 = I.asDiagonal().inverse()*(-w0.cross(I.asDiagonal()*w0));
	Vector3d w1	= w0 + 0.5*dt*Aw1;
	//Second order correction for angular velocity
	Vector3d Aw2 = I.asDiagonal().inverse()*(-w1.cross(I.asDiagonal()*w1));
	W	= w1 + 0.5*dt*Aw2;

	//store the acceleration for next update
	Avb	= Av;
	Awb	= Aw0+Aw1+Aw2;

	//Update the range of warpping box
	MaxX	= (int) (X(0)+R+1);
	MaxY	= (int) (X(1)+R+1);
	MaxZ	= (int) (X(2)+R+1);

	MinX	= (int) (X(0)-R);
	MinY	= (int) (X(1)-R);
	MinZ	= (int) (X(2)-R);
}

inline void DEM_PARTICLE::ZeroForceTorque()
{
	//set zero hydro and contact force and torque after moved.
	Fh 	<< 0, 0, 0;
	Fc 	<< 0, 0, 0;
	Th 	<< 0, 0, 0;
	Tc 	<< 0, 0, 0;
}

class SPHERE:public DEM_PARTICLE
{
public:
	SPHERE(int tag, double r, const Vector3d& x, double rho);
};

inline SPHERE::SPHERE(int tag, double r, const Vector3d& x, double rho):DEM_PARTICLE(tag, x, rho)
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

	MaxX	= (int) (X(0)+R+1);
	MaxY	= (int) (X(1)+R+1);
	MaxZ	= (int) (X(2)+R+1);

	MinX	= (int) (X(0)-R);
	MinY	= (int) (X(1)-R);
	MinZ	= (int) (X(2)-R);
}

class DISK2D:public DEM_PARTICLE
{
public:
	DISK2D(int tag, double r, const Vector3d& x, double rho);
};

inline DISK2D::DISK2D(int tag, double r, const Vector3d& x, double rho):DEM_PARTICLE(tag, x, rho)
{
    Type= 0;
	R	= r;
	M	= Rho*M_PI*R*R;
	I(0)	= 0.25*M*R*R;
	I(1)	= I(0);
	I(2)	= 2.*I(0);

	P0.resize(1);
	P0[0] << 1., 0., 0.;
	P.resize(1);
	P[0] << 1., 0., 0.;

	MaxX	= (int) (X(0)+R+1);
	MaxY	= (int) (X(1)+R+1);
	MaxZ	= 0;

	MinX	= (int) (X(0)-R);
	MinY	= (int) (X(1)-R);
	MinZ	= 0;
}