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

double Sign(double x)
{
	return (x<0)? -1:1;
}

// ======================================================================
// 1D linear shape function
double ShapeL(double x, double xc, double lx)
{
	double nx = 0.;
	double n = abs(x-xc)/lx;
	if (n<1.)	nx = 1.-n;
	return nx;
}
// 1D derivative of linear shape function
// Not sure the derivative at zero point.
double DShapeL(double x, double xc, double lx)
{
	double d = 0;
	double n = x-xc;
	if (abs(n)<lx)	d = -Sign(n)/lx;
	// if 		(n>=0)	d = -1./lx;
	// else if (n<0)	d =  1./lx;
	return d;
}

void LS1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double nx = ShapeL(x(0), xc(0), l(0));
	double dx = DShapeL(x(0), xc(0), l(0));
	n = nx;
	gn(0) = dx;
	gn(1) = 0.;
	gn(2) = 0.;
}

void LS2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double nx = ShapeL(x(0), xc(0), l(0));
	double ny = ShapeL(x(1), xc(1), l(1));
	double dx = DShapeL(x(0), xc(0), l(0));
	double dy = DShapeL(x(1), xc(1), l(1));
	n = nx*ny;
	gn(0) = dx*ny;
	gn(1) = nx*dy;
	gn(2) = 0.;
}

void LS3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double nx = ShapeL(x(0), xc(0), l(0));
	double ny = ShapeL(x(1), xc(1), l(1));
	double nz = ShapeL(x(2), xc(2), l(2));
	double dx = DShapeL(x(0), xc(0), l(0));
	double dy = DShapeL(x(1), xc(1), l(1));
	double dz = DShapeL(x(2), xc(2), l(2));
	n = nx*ny*nz;
	gn(0) = dx*ny*nz;
	gn(1) = nx*dy*nz;
	gn(2) = nx*ny*dz;
}

// 1D linear shape function
double ShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeL(x(0),xc(0),l(0));
}

// 1D gradient of linear shape function
Vector3d GradShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = DShapeL(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D linear shape function
double ShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1));
}

// 2D gradient of linear shape function
Vector3d GradShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeL(x(1), xc(1), l(1))*DShapeL(x(0), xc(0), l(0));
	grad(1) = ShapeL(x(0), xc(0), l(0))*DShapeL(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3d linear shape function
double ShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1)) * ShapeL(x(2),xc(2),l(2));
}

// 3D gradient of linear shape function
Vector3d GradShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeL(x(1), xc(1), l(1))*ShapeL(x(2), xc(2), l(2))*DShapeL(x(0), xc(0), l(0));
	grad(1) = ShapeL(x(0), xc(0), l(0))*ShapeL(x(2), xc(2), l(2))*DShapeL(x(1), xc(1), l(1));
	grad(2) = ShapeL(x(0), xc(0), l(0))*ShapeL(x(1), xc(1), l(1))*DShapeL(x(2), xc(2), l(2));

	return grad;
}

// ======================================================================
// 1D quadratic B-spline shape function
double ShapeQ(double x, double xc, double lx)
{
	double nx = 0.;
	double n = abs(x - xc)/lx;
	if (n<1.5)
	{
		if (n> 0.5)		nx = 0.5*(1.5-n)*(1.5-n);
		else			nx = 0.75 - n*n;
	}
	return nx;
}
// Derivative of 1D quadratic B-spline shape function
double DShapeQ(double x, double xc, double lx)
{
	double d = 0;
	double s = Sign(x-xc);
	double n = abs(x - xc)/lx;
	if (n<1.5)
	{
		if (n> 0.5)		d = (n - 1.5)*s;
		else			d = -2.*n*s;
	}
	return d;
}
// 1D quadratic B-spline shape function
double ShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeQ(x(0),xc(0),l(0));
}
// 1D gradient of quadratic B-spline shape function
Vector3d GradShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = DShapeQ(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D quadratic B-spline shape function
double ShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1));
}
// 2D gradient of quadratic B-spline shape function
Vector3d GradShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeQ(x(1), xc(1), l(1))*DShapeQ(x(0), xc(0), l(0));
	grad(1) = ShapeQ(x(0), xc(0), l(0))*DShapeQ(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3D quadratic B-spline shape function
double ShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1)) * ShapeQ(x(2),xc(2),l(2));
}
// 3D gradient of quadratic B-spline shape function
Vector3d GradShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeQ(x(1), xc(1), l(1))*ShapeQ(x(2), xc(2), l(2))*DShapeQ(x(0), xc(0), l(0));
	grad(1) = ShapeQ(x(0), xc(0), l(0))*ShapeQ(x(2), xc(2), l(2))*DShapeQ(x(1), xc(1), l(1));
	grad(2) = ShapeQ(x(0), xc(0), l(0))*ShapeQ(x(1), xc(1), l(1))*DShapeQ(x(2), xc(2), l(2));

	return grad;
}
// ======================================================================
// 1D cubic B-spline shape function
double ShapeC(double x, double xc, double lx)
{
	double nx = 0.;
	double n = abs(x - xc)/lx;
	if (n<2.)
	{
		if (n> 1.)		nx = pow(2.-n, 3.)/6.;
		else			nx = 0.5*n*n*n - n*n + 0.66666666666666666;
	}
	return nx;
}
// Derivative of 1D cubic B-spline shape function
double DShapeC(double x, double xc, double lx)
{
	double d = 0;
	double s = Sign(x-xc);
	double n = abs(x - xc)/lx;
	if (n<2.)
	{
		if (n> 1.)		d = -0.5*pow(2.-n, 2.)*s;
		else			d = (1.5*n*n -2.*n)*s;
	}
	return d;
}
// 1D cubic B-spline shape function
double ShapeC1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeC(x(0),xc(0),l(0));
}
// 1D gradient of cubic B-spline shape function
Vector3d GradShapeC1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = DShapeC(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D cubic B-spline shape function
double ShapeC2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeC(x(0),xc(0),l(0)) * ShapeC(x(1),xc(1),l(1));
}
// 2D gradient of cubic B-spline shape function
Vector3d GradShapeC2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeC(x(1), xc(1), l(1))*DShapeC(x(0), xc(0), l(0));
	grad(1) = ShapeC(x(0), xc(0), l(0))*DShapeC(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3D cubic B-spline shape function
double ShapeC3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	return ShapeC(x(0),xc(0),l(0)) * ShapeC(x(1),xc(1),l(1)) * ShapeC(x(2),xc(2),l(2));
}
// 3D gradient of cubic B-spline shape function
Vector3d GradShapeC3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
{
	Vector3d grad;

	grad(0) = ShapeC(x(1), xc(1), l(1))*ShapeC(x(2), xc(2), l(2))*DShapeC(x(0), xc(0), l(0));
	grad(1) = ShapeC(x(0), xc(0), l(0))*ShapeC(x(2), xc(2), l(2))*DShapeC(x(1), xc(1), l(1));
	grad(2) = ShapeC(x(0), xc(0), l(0))*ShapeC(x(1), xc(1), l(1))*DShapeC(x(2), xc(2), l(2));

	return grad;
}
// ======================================================================
// 1D GIMP shape function
/*Only considering uGIMP right now 14 Jun 2018*/
double ShapeGIMP(double x, double xc, double lx, double lpx)
{
	double d = 0;
	double n = x-xc;
	if (n<(lx+lpx))
	{
		if 		(n>(lx-lpx) )	d = pow(lx+lpx-n,2)/(4.*lx*lpx);
		else if (n>lpx      )	d = 1.-n/lx;
		else if (n>(-lpx)   )	d = 1.-(n*n/lpx+lpx)/(2.*lx);
		else if (n>(-lx+lpx))	d = 1.+n/lx;
		else if (n>(-lx-lpx))	d = pow(lx+lpx+n,2)/(4.*lx*lpx);
	}
	return d;
}
// Derivative of 1D GIMP shape function
double DShapeGIMP(double x, double xc, double lx, double lpx)
{
	double d = 0;
	double n = x-xc;
	if (n<(lx+lpx))
	{
		if 		(n>(lx-lpx) )	d = -(lx+lpx-n)/(2.*lx*lpx);
		else if (n>lpx      )	d = -1./lx;
		else if (n>(-lpx)   )	d = -n/(lx*lpx);
		else if (n>(-lx+lpx))	d = 1./lx;
		else if (n>(-lx-lpx))	d = (lx+lpx+n)/(2.*lx*lpx);
	}
	return d;
}
void GIMP(double x, double xc, double lx, double lpx, double& n, double& gn)
{
	n = 0.;
	gn = 0.;
	double d = x-xc;
	double a = lx+lpx;
	double b = lx*lpx;
	if (d<a)
	{
		if (d>(lx-lpx))
		{
			n 	= (a-d)*(a-d)/(4.*b);
			gn 	= (d-a)/(2.*b);
		}
		else if (d>lpx)
		{
			n 	= 1.-d/lx;
			gn 	= -1./lx;
		}
		else if (d>(-lpx))
		{
			n 	= 1.-(d*d+lpx*lpx)/(2.*b);
			gn 	= -d/b;
		}
		else if (d>(-lx+lpx))
		{
			n 	= 1.+d/lx;
			gn 	= 1./lx;
		}
		else if (d>(-a))
		{
			n 	= (a+d)*(a+d)/(4.*b);
			gn 	= (a+d)/(2.*b);
		}
	}
}

void GIMP1D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double n0 = 0.;
	double gn0 = 0.;

	GIMP(x(0), xc(0), l(0), lp(0), n0, gn0);

	n = n0;

	gn(0) = gn0;
	gn(1) = 0.;
	gn(2) = 0.;
}

void GIMP2D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double n0 = 0.;
	double n1 = 0.;
	double gn0 = 0.;
	double gn1 = 0.;

	GIMP(x(0), xc(0), l(0), lp(0), n0, gn0);
	GIMP(x(1), xc(1), l(1), lp(1), n1, gn1);

	n = n0*n1;

	gn(0) = gn0*n1;
	gn(1) = gn1*n0;
	gn(2) = 0.;
}

void GIMP3D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
{
	double n0 = 0.;
	double n1 = 0.;
	double n2 = 0.;
	double gn0 = 0.;
	double gn1 = 0.;
	double gn2 = 0.;

	GIMP(x(0), xc(0), l(0), lp(0), n0, gn0); 
	GIMP(x(1), xc(1), l(1), lp(1), n1, gn1);
	GIMP(x(2), xc(2), l(2), lp(2), n2, gn2);

	n = n0*n1*n2;

	gn(0) = gn0*n1*n2;
	gn(1) = gn1*n0*n2;
	gn(2) = gn2*n0*n1;
}
// ======================================================================
// MLS
// polynomials
VectorXd PQ2D(Vector3d& x, Vector3d& xc)
{
	VectorXd p(6);
	// 1, x, y, x^2, xy, y^2
	p(0) = 1.;
	p(1) = x(0)-xc(0);
	p(2) = x(1)-xc(1);
	p(3) = p(1)*p(1);
	p(4) = p(1)*p(2);
	p(5) = p(2)*p(2);

	// p(1) = x(0);
	// p(2) = x(1);
	// p(3) = x(0)*x(0);
	// p(4) = x(0)*x(1);
	// p(5) = x(1)*x(1);

	return p;
}
// quadratic spline weight function
double WeightQ(double r, double h)
{
	double w = 0.;
	if (r<0.5*h)		w = 0.75-r*r/(h*h);
	else if (r<1.5*h)	w = 0.5*pow(1.5-r/h,2);
	return w;
}
// 3D quadratic spline weight function
double WeightQ2D(Vector3d& x, Vector3d& xc, Vector3d& l)
{
	double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1));
	return w;
}
// 3D quadratic spline weight function
double WeightQ3D(Vector3d& x, Vector3d& xc, Vector3d& l)
{
	double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1))*WeightQ(abs(x(2)-xc(2)), l(2));
	return w;
}
// cubic spline weight function (need mofided for h)
double WeightC(double r, double h)
{
	double w = 0.;
	if (r<0.5)		w = 0.66666666666+4.*r*r*(r-1.);
	else if (r<1.)	w = 0.66666666666+4.*r*(-1.+r-r*r/3.);
	return w;
}
// 3D cubic spline weight function
double WeightC2D(Vector3d& x, Vector3d& xc, Vector3d& l)
{
	double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1));
	return w;
}
// 3D cubic spline weight function
double WeightC3D(Vector3d& x, Vector3d& xc, Vector3d& l)
{
	double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1))*WeightQ(abs(x(2)-xc(2)), l(2));
	return w;
}

VectorXd MLS(vector<Vector3d>& Xp, Vector3d& xc, size_t wtype)
{
	size_t n = Xp.size();

	if (n>800)
	{
		cout << "n= " << n << endl;
		abort();
	}

	MatrixXd Pm(n,6);
	MatrixXd Wm(n,n);
	Wm.setZero();
	for (size_t i=0; i<n; ++i)
	{
		VectorXd pi = PQ2D(Xp[i], xc);
		Pm.row(i) = pi.transpose();
		// cout << pi.transpose() << endl;
		// Wm(i,i) = WeightC(abs(Xp[i](0)-xc(0)), 1.)*WeightC(abs(Xp[i](1)-xc(1)), 1.);
		if (wtype==0) Wm(i,i) = WeightQ(abs(Xp[i](0)-xc(0)), 1.)*WeightQ(abs(Xp[i](1)-xc(1)), 1.);
		else if (wtype==1) Wm(i,i) = WeightC(abs(Xp[i](0)-xc(0)), 1.)*WeightC(abs(Xp[i](1)-xc(1)), 1.);
		// if (wtype==0) Wm(i,i) = WeightQ((Xp[i]-xc).norm(), 1.);
		// else if (wtype==1) Wm(i,i) = WeightC((Xp[i]-xc).norm(), 1.);
	}
	// cout << "n= " << n << endl;
	// Eq. 10
	MatrixXd Mm = Pm.transpose()*(Wm*Pm);
	// additional constraints
	// if (n>2)
	// {
		Mm(5,5) += 1.0e-2;
		Mm(4,4) += 1.0e-2;
		Mm(3,3) += 1.0e-2;
	// }
	// Eq. 11
	VectorXd Pc = PQ2D(xc,xc);
	// cout << "start phi" << endl;
	VectorXd phi = Pc.transpose()*Mm.inverse()*Pm.transpose()*Wm;
	// cout << "finish phi" << endl;
	// VectorXd phi;
	return phi;
}

VectorXd MLS1(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype)
{
	size_t n = Xp.size();
	MatrixXd Pm(n,6);
	MatrixXd Wm(n,n);
	Wm.setZero();
	for (size_t i=0; i<n; ++i)
	{
		VectorXd pi = PQ2D(Xp[i], xc);
		Pm.row(i) = pi.transpose();
		// cout << pi.transpose() << endl;
		if (wtype==0) Wm(i,i) = WeightQ(abs(Xp[i](0)-xc(0)), 1.)*WeightQ(abs(Xp[i](1)-xc(1)), 1.);
		else if (wtype==1) Wm(i,i) = WeightC(abs(Xp[i](0)-xc(0)), 1.)*WeightC(abs(Xp[i](1)-xc(1)), 1.);
		cout << "Wm(i,i)= " << Wm(i,i) << endl;
		cout << "WeightQ(abs(Xp[i](0)-xc(0)), 1.)= " << WeightQ(abs(Xp[i](0)-xc(0)), 1.) << endl;
		cout << "WeightQ(abs(Xp[i](1)-xc(1)), 1.)= " << WeightQ(abs(Xp[i](1)-xc(1)), 1.) << endl;
		cout << "abs(Xp[i](1)-xc(1))= " << abs(Xp[i](1)-xc(1)) << endl;
		cout << "Xp[i](1)= " << Xp[i](1) << endl;
		cout << "xc(1)= " << xc(1) << endl;
		cout << "Xp[i](1)-xc(1)= " << Xp[i](1)-xc(1) << endl;

	}
	// Eq. 10
	MatrixXd Mm = Pm.transpose()*Wm*Pm;
	// additional constraints
	Mm(5,5) += vis;
	Mm(4,4) += vis;
	Mm(3,3) += vis;
	// Eq. 11
	VectorXd Pc = PQ2D(xc,xc);
	VectorXd phi = Pc.transpose()*Mm.inverse()*Pm.transpose()*Wm;
	return phi;
}

VectorXd MLS2(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype)
{
	size_t n = Xp.size();
	MatrixXd Pm(n,6);
	MatrixXd Wm(n,n);
	Wm.setZero();
	for (size_t i=0; i<n; ++i)
	{
		VectorXd pi = PQ2D(Xp[i], xc);
		Pm.row(i) = pi.transpose();

		if (wtype==0) Wm(i,i) = WeightQ(abs(Xp[i](0)-xc(0)), 1.)*WeightQ(abs(Xp[i](1)-xc(1)), 1.);
		else if (wtype==1) Wm(i,i) = WeightC(abs(Xp[i](0)-xc(0)), 1.)*WeightC(abs(Xp[i](1)-xc(1)), 1.);
		cout << "Wm(i,i)= " << Wm(i,i) << endl;
		cout << "abs(Xp[i](0)-xc(0))= " << abs(Xp[i](0)-xc(0)) << endl;
		cout << "abs(Xp[i](1)-xc(1))= " << abs(Xp[i](1)-xc(1)) << endl;
	}
	// Eq. 10
	MatrixXd Mm = Pm.transpose()*Wm*Pm;
	// additional constraints
	Mm(5,5) += vis;
	Mm(4,4) += vis;
	Mm(3,3) += vis;
	// Eq. 11
	Vector3d x0 (0.,0.,0.);
	VectorXd Pc = PQ2D(xc,xc);
	VectorXd phi = Pc.transpose()*Mm.inverse()*Pm.transpose()*Wm;
	return phi;
}