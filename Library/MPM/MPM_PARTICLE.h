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

class MPM_PARTICLE
{
public:
	MPM_PARTICLE();
	MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson);
	void Elastic(Matrix3d de);
	void Newtonian(Matrix3d de);
	void MohrCoulomb(Matrix3d de);
	void EOSMorris(double C);
	void EOSMonaghan(double C);
	// void DruckerPrager(Matrix3d de);

    int 						Type;                       // Type of particle, for -1 is freely moved particles or -2 is boundary particles.
	int 						ID; 				    	// Index of particle in the list 
	int 						Tag;				    	// Tag of particle

	double 						M;				            // Mass
	double 						Vol;						// Volume
	double 						Vol0;						// Init Volume
	double 						Arc;						// Arc length for boundary nodes

	double 						Mu;							// Shear modulus (Lame's first parameter)
	double						La;							// Lame's second parameter

	double 						Young;						// Young's modus
	double						Poisson;					// Possion ratio
	double 						C;							// Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double 						Phi;						// Angle of internal friction
	double 						Psi;						// Angle of dilatation

	double 						P;							// Pressure of fluid

	Vector3d					PSize0;						// Vector of half length of particle domain at init
	Vector3d					PSize;						// Vector of half length of particle domain

	Vector3d 					X;				            // Position
	Vector3d 					X0;				            // Init position
	Vector3d 					DeltaX;				        // Increasement of position
	Vector3d					V;							// Velocity
	Vector3d					Vf;							// Fixed velocity
	Vector3d					B;							// Body force acc
	Vector3d					Fh0;						// Hydro force
	Vector3d					Fh;							// Hydro force
	Vector3d					Nor;						// Normal direction (only non-zero for boundary particles)

	Matrix3d					S;							// Stress
	Matrix3d					L;							// Velocity gradient tensor
	Matrix3d					F;							// Derformation gradient tensor
	Matrix3d					Dp;							// Elastic tensor in principal stress space
	Matrix3d					Dpi;						// Inverse of Dp

	bool						FixV;						// Whether the velocity is fixed
	bool						Removed;					// whether this particle is removed

	vector<int>					Lnei;						// List of neighor nodes indexs, used to calculate arc lengh for FSI problems
	vector<Vector3i>			Ln;							// List of node indexs
	vector<double>				LnN;						// List of shape functions
	vector<Vector3d>			LnGN;						// List of gradient of shape functions
};

inline MPM_PARTICLE::MPM_PARTICLE()
{
    Type	= -1;
	ID		= 0;
	Tag		= 0;
	M 		= 0.;
	X 		= Vector3d::Zero();
	X0 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fh0 	= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	S 		= Matrix3d::Zero();
	F 		= Matrix3d::Identity();

	FixV	= false;
	Removed	= false;
}

inline MPM_PARTICLE::MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson)
{
    Type	= type;
	ID		= 0;
	Tag		= 0;
	M 		= m;
	X 		= x;
	X0 		= x;
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	S 		= Matrix3d::Zero();
	F 		= Matrix3d::Identity();

	FixV	= false;
	Removed	= false;

	Young 	= young;
	Poisson = poisson;

	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);

	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;

	Dpi = Dp.inverse();
}

// Elastic model
inline void MPM_PARTICLE::Elastic(Matrix3d de)
{
	S += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();
}

// Elastic model
inline void MPM_PARTICLE::Newtonian(Matrix3d de)
{
	S += 2.*Mu*de - (0.666666666666*Mu-P)*de.trace()*Matrix3d::Identity();
}

void MPM_PARTICLE::EOSMorris(double C)
{
	P = C*C*M/Vol;
}

void MPM_PARTICLE::EOSMonaghan(double C)
{
	P = C*C*M/Vol0/7.*(pow(Vol0/Vol,7.)-1.);
}

// Mohr-Coulomb model
// Based on "An efficient return algorithm for non-associated plasticity with linear yield criteria in principal stress space"
inline void MPM_PARTICLE::MohrCoulomb(Matrix3d de)
{
	// Apply elastic model first
	Elastic(de);
	SelfAdjointEigenSolver<Matrix3d> eigensolver(S);

	double s1 = eigensolver.eigenvalues()(2);
	double s2 = eigensolver.eigenvalues()(1);
	double s3 = eigensolver.eigenvalues()(0);

	Vector3d sb (s1, s2, s3);

	double f = (s1-s3) +(s1+s3)*sin(Phi) -2.*C*cos(Phi);		// Eq.28

	if (f>0.)
	{
		Matrix3d V;

		V.col(0) = eigensolver.eigenvectors().col(2);
		V.col(1) = eigensolver.eigenvectors().col(1);
		V.col(2) = eigensolver.eigenvectors().col(0);

		Vector3d sc;

		double k = (1.+sin(Phi)) / (1.-sin(Phi));				// Eq.32
		Vector3d a1 (k, 0., -1.);								// Eq.32
		double m = (1.+sin(Psi)) / (1.-sin(Psi));				// Eq.33
		Vector3d b1 (m, 0., -1.);								// Eq.33

		Vector3d rp = Dp*b1 / (b1.transpose()*Dp*a1);			// Eq.27b

		Vector3d sa (1., 1., 1.);
		sa *= 2.*C*sqrt(k)/(k-1.);								// Eq.34

		Vector3d rl1 (1., 1., k);								// Eq.40
		Vector3d rl2 (1., k , k);								// Eq.40

		Vector3d rgl1 (1., 1., m);								// Eq.41
		Vector3d rgl2 (1., m , m);								// Eq.41

		double t1 = rgl1.transpose()*Dpi*(sb-sa); 				// Eq.39
		double t2 = rgl2.transpose()*Dpi*(sb-sa);				// Eq.39

		double t1f = rgl1.transpose()*Dpi*rl1;
		double t2f = rgl2.transpose()*Dpi*rl2;

		t1 /= t1f;
		t2 /= t2f;

		double p12 = rp.cross(rl1).dot(sb-sa);					// Eq.45
		double p13 = rp.cross(rl2).dot(sb-sa);					// Eq.46

		// return to apex
		if (t1>0. && t2>0.)
		{
			sc = sa;											// Eq.42
		}
		// return to plane f=0
		else if (p12>=0. && p13<=0.)
		{
			Vector3d dsp = f*rp;								// Eq.27a
			sc = sb - dsp;										// Eq.6
		}
		// return to l1
		else if (p12<0. && p13<0.)
		{
			sc = t1*rl1 + sa;									// Eq.40
		}
		// return to l2
		else if (p12>0. && p13>0.)
		{
			sc = t2*rl2 + sa;									// Eq.40
		}

		Matrix3d sp = Matrix3d::Zero();
		sp(0,0) = sc(0);
		sp(1,1) = sc(1);
		sp(2,2) = sc(2);

		S = V * sp * V.inverse();
	}
}

// void MPM_PARTICLE::DruckerPrager(Matrix3d de)
// {
// 	auto yieldFunc = [](Matrix3d s, double c)
// 	{
// 		double p = s.trace()/3.;
// 		Matrix3d ss = s - p*Matrix3d::Identity();
// 		double j2 = 0.5*(ss.array()*ss.array()).sum();
// 		return sqrt(j2)+a*p-b*c;
// 	};
// 	// Apply elastic model first
// 	Elastic(de);	
// }

// Matsuoka-Nakai model
/*void MPM_PARTICLE::MNModel()
{
	double sinPhi2 = sin(Phi)*sin(Phi);
	double kf = (sinPhi2-9.)/(sinPhi2-1.);

	// Matrix3d Sb = S - C*Matrix3d::Identity();

	double I1 = S.trace();
	double I2 = 0.5*(I1*I1 - (S*S).trace());
	double I3 = S.determinant();

	double f = 6.*(I3*kf -I1*I2) + C*(12.*I1*I1 + 18.*I2 - 6.*I2*kf) + C*C*(6.*I1*kf - 54.*I1) + C*C*C*(54. - 6.*kf);
}*/