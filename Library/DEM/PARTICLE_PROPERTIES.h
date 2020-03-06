// This part of code is used to calculate mass center, volume and inertia tensor of CONVEX polyhedron
// Mainly based on "fast and accurate computation of polyhedral mass properties"
// The original C code can be found at: https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html
// Note: the original code have mirror error (last three lines)

void CalProjectionIntegrals(VectorXi& f, vector<Vector3d>& P, Vector3i& d, VectorXd& pi)
{
	pi.resize(10);
	pi.setZero();
	double a0, a1, da;
	double b0, b1, db;
	double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
	double a1_2, a1_3, b1_2, b1_3;
	double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
	double Cab, Kab, Caab, Kaab, Cabb, Kabb;

	for (int i=0; i<f.size(); ++i)
	{
		a0 = P[f(i)](d(0));
		b0 = P[f(i)](d(1));
		a1 = P[f((i+1)%f.size())](d(0));
		b1 = P[f((i+1)%f.size())](d(1));
	    da = a1 - a0;
	    db = b1 - b0;
	    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
	    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
	    a1_2 = a1 * a1; a1_3 = a1_2 * a1; 
	    b1_2 = b1 * b1; b1_3 = b1_2 * b1;

	    C1 = a1 + a0;
	    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
	    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
	    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
	    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
	    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
	    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

	    pi(0) += db*C1;
	    pi(1) += db*Ca;
	    pi(2) += db*Caa;
	    pi(3) += db*Caaa;
	    pi(4) += da*Cb;
	    pi(5) += da*Cbb;
	    pi(6) += da*Cbbb;
	    pi(7) += db*(b1*Cab + b0*Kab);
	    pi(8) += db*(b1*Caab + b0*Kaab);
	    pi(9) += da*(a1*Cabb + a0*Kabb);
	}
	pi(0) /= 2.0;
	pi(1) /= 6.0;
	pi(2) /= 12.0;
	pi(3) /= 20.0;
	pi(4) /= -6.0;
	pi(5) /= -12.0;
	pi(6) /= -20.0;
	pi(7) /= 24.0;
	pi(8) /= 60.0;
	pi(9) /= -60.0;
}

void CalFaceIntegrals(VectorXi& f, vector<Vector3d>& P, Vector3i& d, VectorXd& fi)
{
	fi.resize(12);

	Vector3d n = ((P[f(2)]-P[f(1)]).cross(P[f(0)]-P[f(1)])).normalized();
	double w = -n.dot(P[f(0)]);

	double k1, k2, k3, k4;
	k1 = 1 / n(d(2)); k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

	VectorXd pi;
	CalProjectionIntegrals(f, P, d, pi);

	Vector3d sqrn (n(d(0))*n(d(0)), n(d(1))*n(d(1)), 0.);
	Vector3d cuben (n(d(0))*n(d(0))*n(d(0)), n(d(1))*n(d(1))*n(d(1)), 0.);

	fi(0) = k1 * pi(1);
	fi(1) = k1 * pi(4);
	fi(2) = -k2 * (n(d(0))*pi(1) + n(d(1))*pi(4) + w*pi(0));

	fi(3) = k1 * pi(2);
	fi(4) = k1 * pi(5);
	fi(5) = k3 * (sqrn(0)*pi(2) + 2*n(d(0))*n(d(1))*pi(7) + sqrn(1)*pi(5)
	 + w*(2*(n(d(0))*pi(1) + n(d(1))*pi(4)) + w*pi(0)));

	fi(6) = k1 * pi(3);
	fi(7) = k1 * pi(6);
	fi(8) = -k4 * (cuben(0)*pi(3) + 3*sqrn(0)*n(d(1))*pi(8) 
	   + 3*n(d(0))*sqrn(1)*pi(9) + cuben(1)*pi(6)
	   + 3*w*(sqrn(0)*pi(2) + 2*n(d(0))*n(d(1))*pi(7) + sqrn(1)*pi(5))
	   + w*w*(3*(n(d(0))*pi(1) + n(d(1))*pi(4)) + w*pi(0)));

	fi(9) = k1 * pi(8);
	fi(10) = -k2 * (n(d(0))*pi(9) + n(d(1))*pi(6) + w*pi(5));
	fi(11) = k3 * (sqrn(0)*pi(3) + 2*n(d(0))*n(d(1))*pi(8) + sqrn(1)*pi(9)
	 + w*(2*(n(d(0))*pi(2) + n(d(1))*pi(7)) + w*pi(1)));
}

void CalVolumeIntegrals(vector<VectorXi>& F, vector<Vector3d>& P, double rho, double& vol, double& m, Vector3d& xc, Matrix3d& Ip)
{
	double t0 = 0.;
	Matrix3d ti;
	ti.setZero();

	double nx, ny, nz; 
	for (size_t i=0; i<F.size(); ++i)
	{
		VectorXi f = F[i];
		Vector3d n = ((P[f(2)]-P[f(1)]).cross(P[f(0)]-P[f(1)])).normalized();
		nx = abs(n(0));
		ny = abs(n(1));
		nz = abs(n(2));
		Vector3i d (0,0,0);
	    if (nx > ny && nx > nz) d(2) = 0;
	    else d(2) = (ny > nz) ? 1 : 2;
	    d(0) = (d(2)+1)%3;
	    d(1) = (d(0)+1)%3;

	    VectorXd fi;
	    CalFaceIntegrals(f, P, d, fi);

	    t0 += n(0) * ((d(0) == 0) ? fi(0) : ((d(1) == 0) ? fi(1) : fi(2)));

	    ti(d(0),0) += n(d(0)) * fi(3);
	    ti(d(1),0) += n(d(1)) * fi(4);
	    ti(d(2),0) += n(d(2)) * fi(5);
	    ti(d(0),1) += n(d(0)) * fi(6);
	    ti(d(1),1) += n(d(1)) * fi(7);
	    ti(d(2),1) += n(d(2)) * fi(8);
	    ti(d(0),2) += n(d(0)) * fi(9);
	    ti(d(1),2) += n(d(1)) * fi(10);
	    ti(d(2),2) += n(d(2)) * fi(11);
	}
	ti.col(0) /= 2.;
	ti.col(1) /= 3.;
	ti.col(2) /= 2.;

	vol = t0;
	m = rho*vol;

	xc = ti.col(0)/t0;

	Matrix3d I;

	I(0,0) = rho*(ti(1,1) + ti(2,1));
	I(1,1) = rho*(ti(2,1) + ti(0,1));
	I(2,2) = rho*(ti(0,1) + ti(1,1));
	I(0,1) = I(1,0) = rho*ti(0,2);
	I(2,1) = I(1,2) = rho*ti(1,2);
	I(0,2) = I(2,0) = rho*ti(2,2);

	I(0,0) -= m*(xc(1)*xc(1) + xc(2)*xc(2));
	I(1,1) -= m*(xc(2)*xc(2) + xc(0)*xc(0));
	I(2,2) -= m*(xc(0)*xc(0) + xc(1)*xc(1));
	I(0,1) = I(1,0) -= m*xc(0)*xc(1);
	I(2,1) = I(1,2) -= m*xc(1)*xc(2);
	I(0,2) = I(2,0) -= m*xc(0)*xc(2);

	Ip = I;
}