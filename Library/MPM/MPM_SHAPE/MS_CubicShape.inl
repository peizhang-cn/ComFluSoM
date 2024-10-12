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

namespace MPM_ShapeFunction
{
	// ======================================================================
	// 1D cubic B-spline shape function
	inline double ShapeCubicBSpline1(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = abs(x - xc)/lx;
		if (n<2.)
		{
			double n3 = n*n*n;
			if (n> 1.)		nx = -n3/6. + n*n -2.*n + 4./3.;
			else			nx = n3/6. - n + 1.;
		}
		return nx;
	}

	inline double ShapeCubicBSpline2(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<2.)
		{
			double n2 = n*n;
			double n3 = n2*n;
			if (n> 1.)			nx = -n3/6. + n2 -2.*n + 4./3.;
			else if (n>0.)		nx = 0.5*n3 - n2 + 2./3.;
			else if (n>-1.)		nx = -n3/3. - n2 + 2./3.;
		}
		return nx;
	}	

	inline double ShapeCubicBSpline3(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = abs(x - xc)/lx;
		if (n<2.)
		{
			if (n> 1.)		nx = pow(2.-n, 3.)/6.;
			else			nx = 0.5*n*n*n - n*n + 2./3.;
		}
		return nx;
	}

	inline double ShapeCubicBSpline4(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<1.)
		{
			double n2 = n*n;
			double n3 = n2*n;
			if (n> 0.)			nx = n3/3. - n2 + 2./3.;
			else if (n>-1.)		nx = -0.5*n3 - n2 + 2./3.;
			else if (n>-2.)		nx = n3/6. + n2 + 2.*n + 4./3.;
		}
		return nx;
	}

	// Derivative of 1D cubic B-spline shape function
	inline double DShapeCubicBSpline1(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<2.)
		{
			double n2 = n*n;
			if (n>1.)		nx = -0.5*n2 + 2.*n - 2.;
			else if (n>0.)	nx = 0.5*n2 -1.;
			else if (n>-1.)	nx = -0.5*n2 + 1.;
			else if (n>-2.)	nx = 0.5*n2 + 2.*n + 2.;
		}
		return nx;
	}

	inline double DShapeCubicBSpline2(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<2.)
		{
			double n2 = n*n;
			if (n>1.)		nx = -0.5*n2 + 2.*n - 2.;
			else if (n>0.)	nx = 1.5*n2 -2.*n;
			else if (n>-1.)	nx = -n2 -2.*n;
		}
		return nx;
	}

	inline double DShapeCubicBSpline3(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<2.)
		{
			double n2 = n*n;
			if (n>1.)		nx = -0.5*n2 + 2.*n - 2.;
			else if (n>0.)	nx = 1.5*n2 -1.;
			else if (n>-1.)	nx = -1.5*n2 - 1.;
			else if (n>-2.)	nx = 0.5*n2 + 2.*n + 2.;
		}
		return nx;
	}

	inline double DShapeCubicBSpline4(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = (x - xc)/lx;
		if (n<2.)
		{
			double n2 = n*n;
			if (n>1.)		nx = -0.5*n2 + 2.*n - 2.;
			else if (n>0.)	nx = 0.5*n2 -1.;
			else if (n>-1.)	nx = -0.5*n2 + 1.;
			else if (n>-2.)	nx = 0.5*n2 + 2.*n + 2.;
		}
		return nx;
	}

	inline double ShapeCubicBSpline(int type, double x, double xc, double lx)
	{
		double nx = 0.;
		if (type==1)		nx = ShapeCubicBSpline1(x,xc,lx);
		else if (type==2)	nx = ShapeCubicBSpline2(x,xc,lx);
		else if (type==3)	nx = ShapeCubicBSpline3(x,xc,lx);
		else if (type==4)	nx = ShapeCubicBSpline4(x,xc,lx);
		return nx;
	}

	inline double DShapeCubicBSpline(int type, double x, double xc, double lx)
	{
		double nx = 0.;
		if (type==1)		nx = DShapeCubicBSpline1(x,xc,lx);
		else if (type==2)	nx = DShapeCubicBSpline2(x,xc,lx);
		else if (type==3)	nx = DShapeCubicBSpline3(x,xc,lx);
		else if (type==4)	nx = DShapeCubicBSpline4(x,xc,lx);
		return nx;
	}

	// 1D cubic B-spline shape function
	inline double ShapeCubicBSpline1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		int i = (int) type(0);
		return ShapeCubicBSpline(i, x(0),xc(0),l(0));
	}
	// 1D gradient of cubic B-spline shape function
	inline Vector3d GradShapeCubicBSpline1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		Vector3d grad;
		int i = (int) type(0);
		grad(0) = DShapeCubicBSpline(i, x(0), xc(0), l(0));
		grad(1) = 0;
		grad(2) = 0;

		return grad;
	}

	inline void CubicBSpline1D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn)
	{
		int i = (int) type(0);
		double nx = ShapeCubicBSpline(i, x(0),xc(0),l(0));
		Vector3d grad(0.,0.,0.);
		grad(0) = DShapeCubicBSpline(i, x(0), xc(0), l(0));
		n = nx;
		gn = grad;
	}
	// 2D cubic B-spline shape function
	inline double ShapeCubicBSpline2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		int i = (int) type(0);
		int j = (int) type(1);
		return ShapeCubicBSpline(i,x(0),xc(0),l(0)) * ShapeCubicBSpline(j,x(1),xc(1),l(1));
	}
	// 2D gradient of cubic B-spline shape function
	inline Vector3d GradShapeCubicBSpline2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		Vector3d grad;
		int i = (int) type(0);
		int j = (int) type(1);
		grad(0) = ShapeCubicBSpline(j,x(1), xc(1), l(1))*DShapeCubicBSpline(i,x(0), xc(0), l(0));
		grad(1) = ShapeCubicBSpline(i,x(0), xc(0), l(0))*DShapeCubicBSpline(j,x(1), xc(1), l(1));
		grad(2) = 0;

		return grad;
	}
	inline void CubicBSpline2D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn)
	{
		int i = (int) type(0);
		int j = (int) type(1);
		double nx = ShapeCubicBSpline(i, x(0),xc(0),l(0));
		double ny = ShapeCubicBSpline(j, x(1),xc(1),l(1));
		Vector3d grad(0.,0.,0.);
		grad(0) = ny*DShapeCubicBSpline(i,x(0), xc(0), l(0));
		grad(0) = nx*DShapeCubicBSpline(j,x(1), xc(1), l(1));
		n = nx*ny;
		gn = grad;
	}
	// 3D cubic B-spline shape function
	inline double ShapeCubicBSpline3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		int i = (int) type(0);
		int j = (int) type(1);
		int k = (int) type(2);
		return ShapeCubicBSpline(i,x(0),xc(0),l(0)) * ShapeCubicBSpline(j,x(1),xc(1),l(1)) * ShapeCubicBSpline(k,x(2),xc(2),l(2));
	}
	// 3D gradient of cubic B-spline shape function
	Vector3d GradShapeCubicBSpline3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type)
	{
		Vector3d grad;
		int i = (int) type(0);
		int j = (int) type(1);
		int k = (int) type(2);
		double nx = ShapeCubicBSpline(i, x(0), xc(0), l(0));
		double ny = ShapeCubicBSpline(j ,x(1), xc(1), l(1));
		double nz = ShapeCubicBSpline(k, x(2), xc(2), l(2));
		grad(0) = ny*nz*DShapeCubicBSpline(i, x(0), xc(0), l(0));
		grad(1) = nx*nz*DShapeCubicBSpline(j, x(1), xc(1), l(1));
		grad(2) = nx*ny*DShapeCubicBSpline(k, x(2), xc(2), l(2));

		return grad;
	}

	inline void CubicBSpline3D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn)
	{
		int i = (int) type(0);
		int j = (int) type(1);
		int k = (int) type(2);
		double nx = ShapeCubicBSpline(i, x(0), xc(0), l(0));
		double ny = ShapeCubicBSpline(j ,x(1), xc(1), l(1));
		double nz = ShapeCubicBSpline(k, x(2), xc(2), l(2));
		n = nx*ny*nz;
		Vector3d grad(0.,0.,0.);
		grad(0) = ny*nz*DShapeCubicBSpline(i, x(0), xc(0), l(0));
		grad(1) = nx*nz*DShapeCubicBSpline(j, x(1), xc(1), l(1));
		grad(2) = nx*ny*DShapeCubicBSpline(k, x(2), xc(2), l(2));
		gn = grad;
	}
}