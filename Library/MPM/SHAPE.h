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
	if 		(n>=0)	d = -1./lx;
	else if (n<0)	d =  1./lx;
	return d;
}
// 1D linear shape function
double ShapeL1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeL(x(0),xc(0),l(0));
}

// 1D gradient of linear shape function
Vector3d GradShapeL1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = DShapeL(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D linear shape function
double ShapeL2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1));
}

// 2D gradient of linear shape function
Vector3d GradShapeL2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = ShapeL(x(1), xc(1), l(1))*DShapeL(x(0), xc(0), l(0));
	grad(1) = ShapeL(x(0), xc(0), l(0))*DShapeL(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3d linear shape function
double ShapeL3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1)) * ShapeL(x(2),xc(2),l(2));
}

// 3D gradient of linear shape function
Vector3d GradShapeL3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
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
double ShapeQ1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeQ(x(0),xc(0),l(0));
}
// 1D gradient of quadratic B-spline shape function
Vector3d GradShapeQ1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = DShapeQ(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D quadratic B-spline shape function
double ShapeQ2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1));
}
// 2D gradient of quadratic B-spline shape function
Vector3d GradShapeQ2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = ShapeQ(x(1), xc(1), l(1))*DShapeQ(x(0), xc(0), l(0));
	grad(1) = ShapeQ(x(0), xc(0), l(0))*DShapeQ(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3D quadratic B-spline shape function
double ShapeQ3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1)) * ShapeQ(x(2),xc(2),l(2));
}
// 3D gradient of quadratic B-spline shape function
Vector3d GradShapeQ3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
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
double ShapeC1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeC(x(0),xc(0),l(0));
}
// 1D gradient of cubic B-spline shape function
Vector3d GradShapeC1D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = DShapeC(x(0), xc(0), l(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D cubic B-spline shape function
double ShapeC2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeC(x(0),xc(0),l(0)) * ShapeC(x(1),xc(1),l(1));
}
// 2D gradient of cubic B-spline shape function
Vector3d GradShapeC2D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = ShapeC(x(1), xc(1), l(1))*DShapeC(x(0), xc(0), l(0));
	grad(1) = ShapeC(x(0), xc(0), l(0))*DShapeC(x(1), xc(1), l(1));
	grad(2) = 0;

	return grad;
}
// 3D cubic B-spline shape function
double ShapeC3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeC(x(0),xc(0),l(0)) * ShapeC(x(1),xc(1),l(1)) * ShapeC(x(2),xc(2),l(2));
}
// 3D gradient of cubic B-spline shape function
Vector3d GradShapeC3D (Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
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
// 1D GIMP shape function
double ShapeGIMP1D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeGIMP(x(0),xc(0),l(0), lp(0));
}
// 1D gradient of GIMP shape function
Vector3d GradShapeGIMP1D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = DShapeGIMP(x(0), xc(0), l(0), lp(0));
	grad(1) = 0;
	grad(2) = 0;

	return grad;
}
// 2D GIMP shape function
double ShapeGIMP2D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeGIMP(x(0),xc(0),l(0), lp(0)) * ShapeGIMP(x(1),xc(1),l(1), lp(1));
}
// 2D gradient of GIMP shape function
Vector3d GradShapeGIMP2D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = ShapeGIMP(x(1), xc(1), l(1), lp(1))*DShapeGIMP(x(0), xc(0), l(0), lp(0));
	grad(1) = ShapeGIMP(x(0), xc(0), l(0), lp(0))*DShapeGIMP(x(1), xc(1), l(1), lp(1));
	grad(2) = 0;

	return grad;
}
// 3D GIMP shape function
double ShapeGIMP3D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	return ShapeGIMP(x(0), xc(0), l(0), lp(0)) * ShapeGIMP(x(1), xc(1), l(1), lp(1)) * ShapeGIMP(x(2), xc(2), l(2), lp(2));
}
// 3D gradient of GIMP shape function
Vector3d GradShapeGIMP3D(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp)
{
	Vector3d grad;

	grad(0) = ShapeGIMP(x(1), xc(1), l(1), lp(1))*ShapeGIMP(x(2), xc(2), l(2), lp(2))*DShapeGIMP(x(0), xc(0), l(0), lp(0));
	grad(1) = ShapeGIMP(x(0), xc(0), l(0), lp(0))*ShapeGIMP(x(2), xc(2), l(2), lp(2))*DShapeGIMP(x(1), xc(1), l(1), lp(1));
	grad(2) = ShapeGIMP(x(0), xc(0), l(0), lp(0))*ShapeGIMP(x(1), xc(1), l(1), lp(1))*DShapeGIMP(x(2), xc(2), l(2), lp(2));

	return grad;
}