/***********************************************************************
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
	// MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson);
	MPM_PARTICLE(int tag, const Vector3d& x, double m);
	void Elastic(Matrix3d& de);
	void SetElastic(double young, double poisson);
	void SetNewtonian(double miu);
	void SetGranular(double rhop, double d);
	void SetMohrCoulomb(double young, double poisson, double phi, double psi, double c);
	void SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c);
	void SetTensionCutoff(double pmax);
	void Newtonian(Matrix3d& de);
	void Granular(Matrix3d& de);
	void MohrCoulomb(Matrix3d& de);
	void DruckerPrager(Matrix3d& de);
	void EOSMorris(double C);
	void EOSMonaghan(double C);
	// void DruckerPrager(Matrix3d de);

    int 						Type;                       // Type of particle, 0 for elastic 1 for fluid 2 for soil.
	int 						ID; 				    	// Index of particle in the list 
	int 						Tag;				    	// Tag of particle

	double 						M;				            // Mass
	double 						Vol;						// Volume
	double 						Vol0;						// Init Volume
	double 						R;							// Radius for DEMPM
	double 						Arc;						// Arc length for boundary nodes

	double 						Mu;							// Shear modulus (Lame's second parameter) or viscosity for fluid
	double						La;							// Lame's first parameter
	double 						K;							// bulk modulus
	double 						H;							// Liniear hardening modulus
	double 						Young;						// Young's modus
	double						Poisson;					// Possion ratio
	double 						C;							// Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double 						Phi;						// Angle of internal friction
	double 						Psi;						// Angle of dilatation
	double 						A_dp;						// Drucker–Prager parameters
	double 						B_dp;						// Drucker–Prager parameters
	double 						Ad_dp;						// Drucker–Prager parameters
	double 						Pmax;						// Drucker–Prager parameters for tension cutoff (max pressure)
	bool 						TensionCut;					// Drucker–Prager parameters for tension cutoff

	double						Porosity;
	double 						P;							// Pressure of fluid or granular materials
// == parameters for μ(I) rheology ======================================================
	double 						Rhop;						// Density of granular particles (not macro density of granular materials)
	double 						Diameterp;					// Diameter of granular particles
	double 						Mius;
	double 						Miud;
	double 						I0;
// ======================================================================================
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
	Vector3d					Fc;							// Contact force
	Vector3d					Nor;						// Normal direction (only non-zero for boundary particles)

	Matrix3d					Strain;						// Strain
	Matrix3d					StrainP;					// Plastic strain tensor
	Matrix3d					Stress;						// Stress
	Matrix3d					StressSmooth;				// Smoothed Stress for visualization
	Matrix3d					L;							// Velocity gradient tensor
	Matrix3d					Td;							// Derformation gradient tensor
	Matrix3d					Dp;							// Elastic tensor in principal stress space
	Matrix3d					Dpi;						// Inverse of Dp

	bool						FixV;						// Whether the velocity is fixed
	bool						Removed;					// whether this particle is removed

	vector<int>					Lnei;						// List of neighor nodes indexs, used to calculate arc lengh for FSI problems
	vector<size_t>				Lni;						// List of node indexs
	vector<size_t>				Lgi;						// List of gauss point indexs
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
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

// inline MPM_PARTICLE::MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson)
// {
//     Type	= type;
// 	ID		= 0;
// 	Tag		= 0;
// 	M 		= m;
// 	X 		= x;
// 	X0 		= x;
// 	V 		= Vector3d::Zero();
// 	Vf 		= Vector3d::Zero();
// 	B 		= Vector3d::Zero();
// 	Fh 		= Vector3d::Zero();
// 	Nor 	= Vector3d::Zero();

// 	Strain 	= Matrix3d::Zero();
// 	StrainP = Matrix3d::Zero();
// 	Stress 	= Matrix3d::Zero();
// 	F 		= Matrix3d::Identity();

// 	Lni.resize(0);
// 	LnN.resize(0);
// 	LnGN.resize(0);

// 	FixV	= false;
// 	Removed	= false;

// 	Young 	= young;
// 	Poisson = poisson;

// 	Mu 		= 0.5*Young/(1.+Poisson);
// 	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
// 	K 		= La+2./3.*Mu;
// 	H 		= 0.;
// 	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
// 	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;

// 	Dpi = Dp.inverse();

// 	// Matrix3d dpi;
// 	// dpi(0,0) = dpi(1,1) = dpi(2,2) = 1/Young;
// 	// dpi(0,1) = dpi(1,0) = dpi(0,2) = dpi(2,0) = dpi(1,2) = dpi(2,1) = -Poisson/Young;

// 	// cout << Dpi << endl;
// 	// cout << "========" << endl;
// 	// cout << dpi << endl;
// 	// abort();
// }

inline MPM_PARTICLE::MPM_PARTICLE(int tag, const Vector3d& x, double m)
{
    Type	= -1;
	ID		= 0;
	Tag		= tag;
	M 		= m;
	X 		= x;
	X0 		= x;
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	Td 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

inline void MPM_PARTICLE::SetElastic(double young, double poisson)
{
	Type 	= 0;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi = Dp.inverse();
}

inline void MPM_PARTICLE::SetNewtonian(double miu)
{
	Type 	= 1;
	Mu 		= miu;
}

inline void MPM_PARTICLE::SetGranular(double rhop, double d)
{
	Type 	= 5;
	Mius 	= 0.38;
	Miud 	= 0.64;
	I0 		= 0.279;
	Rhop 	= rhop;
	Diameterp = d;
}

inline void MPM_PARTICLE::SetMohrCoulomb(double young, double poisson, double phi, double psi, double c)
{
	Type 	= 2;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 	= Dp.inverse();
	Phi 	= phi;
	Psi 	= psi;
	C 		= c;
}

inline void MPM_PARTICLE::SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c)
{
	Type 	= 3;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 	= Dp.inverse();
	Phi 	= phi;
	Psi 	= psi;
	C 		= c;
	TensionCut = false;
	// inner
	if (dptype==0)
	{
		double bot = sqrt(3)*(3.+sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.+sin(Psi)));
	}
	// outer
	else if (dptype==1)
	{
		double bot = sqrt(3)*(3.-sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.-sin(Psi)));
	}
	// plane strain
	else if (dptype==2)
	{
		double bot = sqrt(9.+12*tan(Phi)*tan(Phi));
		A_dp = 3.*tan(Phi)/bot;
		B_dp = 3./bot;
		Ad_dp = 3.*tan(Psi)/sqrt(9.+12*tan(Psi)*tan(Psi));
	}
}
// only works with Drucker Prager (12.3.2019)
inline void MPM_PARTICLE::SetTensionCutoff(double pmax)
{
	TensionCut 	= true;
	Pmax 		= pmax;
}

// Elastic model
inline void MPM_PARTICLE::Elastic(Matrix3d& de)
{
	Stress += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();
}

// Newtonian fluid model
inline void MPM_PARTICLE::Newtonian(Matrix3d& de)
{
	Stress = 2.*Mu*(de - de.trace()/3.*Matrix3d::Identity()) - P*Matrix3d::Identity();
}

// Newtonian fluid model
inline void MPM_PARTICLE::Granular(Matrix3d& de)
{
	// normal of strain rate
	double dnorm = de.determinant();
	// Inertial number
	double I = 2.*Diameterp*dnorm/sqrt(P/Rhop);
	// Macro fraction coefficient
	double miu = Mius + (Miud-Mius)/(I0/I+1.);
	double mu = 0.5*miu*P/dnorm;
	if (dnorm<1.0e-12)	mu = 0.5*miu*P;
	Stress = 2.*mu*(de - de.trace()/3.*Matrix3d::Identity()) - P*Matrix3d::Identity();
	// Matrix3d s = miu*P*de/*/dnorm*/;
	// // if (dnorm==0.)	s = miu*P*de;
	// Stress = s - P*Matrix3d::Identity();
	
	// /*if (I<0.)*/{
	// cout << "miu= " << miu << endl;
	// cout << "I= " << I << endl;
	// cout << "dnorm= " << dnorm << endl;
	// cout << "mu= " << mu << endl;
	// cout << "P= " << P << endl;
	// if (P<0.)	abort();
	// /*abort();*/
	// }
}

void MPM_PARTICLE::EOSMorris(double Cs)
{
	P = Cs*Cs*M/Vol;
}

void MPM_PARTICLE::EOSMonaghan(double Cs)
{
	P = Cs*Cs*M/Vol0/7.*(pow(Vol0/Vol,7.)-1.);
}

inline void MPM_PARTICLE::MohrCoulomb(Matrix3d& de)
{
	// Apply elastic model first
	Elastic(de);
	SelfAdjointEigenSolver<Matrix3d> eigensolver(Stress);

	double s1 = eigensolver.eigenvalues()(2);
	double s2 = eigensolver.eigenvalues()(1);
	double s3 = eigensolver.eigenvalues()(0);

	Vector3d sb (s1, s2, s3);

	double sin0 = sin(Phi);
	double cos0 = cos(Phi);
	double sin1 = sin(Psi);
	// double cos1 = cos(Psi);

	double f = (s1-s3) +(s1+s3)*sin0 -2.*C*cos0;	

	if (f>1.e-18)
	{
		Matrix3d v0;
		v0.col(0) = eigensolver.eigenvectors().col(2);
		v0.col(1) = eigensolver.eigenvectors().col(1);
		v0.col(2) = eigensolver.eigenvectors().col(0);

		Vector3d sc;

		double sin01 = sin0*sin1;
		double qA0 = (8.*Mu/3.-4.*K)*sin01;
		double qA1 = Mu*(1.+sin0)*(1.+sin1);
		double qA2 = Mu*(1.-sin0)*(1.-sin1);
		double qB0 = 2.*C*cos0;

		double gsl = 0.5*(s1-s2)/(Mu*(1.+sin1));
		double gsr = 0.5*(s2-s3)/(Mu*(1.-sin1));
		double gla = 0.5*(s1+s2-2.*s3)/(Mu*(3.-sin1));
		double gra = 0.5*(2.*s1-s2-s3)/(Mu*(3.+sin1));

		double qsA = qA0-4.*Mu*(1.+sin01);
		double qsB = f;
		
		double qlA = qA0-qA1-2.*qA2;
		double qlB = 0.5*(1.+sin0)*(s1+s2)-(1.-sin0)*s3-qB0;

		double qrA = qA0-2.*qA1-qA2;
		double qrB = (1.+sin0)*s1 - 0.5*(1.-sin0)*(s2+s3) - qB0;

		double qaA = -4.*K*sin01;
		double qaB = 2.*(s1+s2+s3)/3.*sin0 - qB0;

		double minslsr = min(gsl,gsr);
		double maxlara = max(gla,gra);

		if (minslsr>0. && qsA*minslsr+qsB<0.)
		{
			double dl = -qsB/qsA;
			double ds0 = -dl*(2.*K-4.*Mu/3.)*sin1;
			sc(0) = s1+ds0-dl*(2.*Mu*(1.+sin1));
			sc(1) = s2+ds0;
			sc(2) = s3+ds0+dl*(2.*Mu*(1.-sin1));
		}
		else if (gsl>0. && gla>=gsl && qlA*gsl+qlB>=0. && qlA*gla+qlB<=0.)
		{
			// return left edge
			double dl = -qlB/qlA;
			double ds0 = dl*(4.*Mu/3.-2.*K)*sin1;
			sc(0) = sc(1) = 0.5*(s1+s2)+ds0 - dl*Mu*(1.+sin1);
			sc(2) = s3+ds0 + 2.*dl*Mu*(1.-sin1);
			// cout << "left edge" << endl;
			// abort();
		}
		else if (gsr>0. && gra>=gsr && qrA*gsr+qrB>=0. && qrA*gra+qrB<=0.)
		{
			double dl = -qrB/qrA;
			double ds0 = dl*(4.*Mu/3.-2.*K)*sin1;
			sc(0) = s1 + ds0 - 2.*dl*Mu*(1.+sin1);
			sc(1) = sc(2) = 0.5*(s2+s3) + ds0 + dl*Mu*(1.-sin1);
			// cout << "right edge" << endl;
		}
		else if (maxlara>0. && qaA*maxlara+qaB>=-1.e-24)
		{
			sc(0) = sc(1) = sc(2) = C/tan(Phi);
			// cout << "apex" << endl;
		}
		else
		{
			cout << "undefined" << endl;
			cout << s1 << " " << s2 << " " << s3 << endl;
			cout << "minslsr:" << minslsr << endl;
			cout << "qsA*minslsr+qsB: " << qsA*minslsr+qsB << endl;
			cout << "===============" << endl;
			cout << "gsl: " << gsl << endl;
			cout << "gla: " << gla << endl;
			cout << "qlA*gsl+qlB: " << qlA*gsl+qlB << endl;
			cout << "qlA*gla+qlB: " << qlA*gla+qlB << endl;
			cout << "===============" << endl;
			cout << "gsr: " << gsr << endl;
			cout << "gra: " << gra << endl;
			cout << "qrA*gsr+qrB: " << qrA*gsr+qrB << endl;
			cout << "qrA*gra+qrB: " << qrA*gra+qrB << endl;	
			cout << "===============" << endl;		
			cout << "maxlara " << maxlara << endl;
			cout << "qaA*maxlara+qaB: " << qaA*maxlara+qaB << endl;
			cout << "f: " << f << endl;
			abort();
		}

		Matrix3d sp = Matrix3d::Zero();
		sp(0,0) = sc(0);
		sp(1,1) = sc(1);
		sp(2,2) = sc(2);

		Stress = v0 * sp * v0.inverse();
		double fa = (sc(0)-sc(2)) +(sc(0)+sc(2))*sin0 -2.*C*cos0;
		if (abs(fa)>1.0e-12)
		{
			cout << "f before: " << f << endl;
			cout << "f after: " << fa << endl;
			cout << sc.transpose() << endl;
			cout << "ID: " << ID << endl;
			cout << "X: " << X.transpose() << endl;
			cout << "V: " << V.transpose() << endl;
			abort();			
		}
	}
}

void MPM_PARTICLE::DruckerPrager(Matrix3d& de)
{
	// Apply elastic model first
	Elastic(de);
	double p = Stress.trace()/3.;
	Matrix3d ss = Stress - p*Matrix3d::Identity();
	double j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
	double f = j2sqr+A_dp*p-B_dp*C;
	if (f>0.)
	{
		// return to cone
		if (A_dp*(p-j2sqr/Mu*K*Ad_dp)-B_dp*C<0.)
		{
			double dl = (j2sqr+A_dp*p-B_dp*C)/(Mu+A_dp*K*Ad_dp);
			Stress -= dl*(Mu/j2sqr*ss+K*Ad_dp*Matrix3d::Identity());
		}
		// return to apex
		else
		{
			Stress = B_dp*C/A_dp*Matrix3d::Identity();
		}
		p = Stress.trace()/3.;
		ss = Stress - p*Matrix3d::Identity();
		j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
		double fb = f;
		f = j2sqr+A_dp*p-B_dp*C;
		if (f>1.0e-8)
		{
			cout << "f before: " << fb << endl;
			cout << "f after: " << f << endl;
			abort();
		}
	}
}