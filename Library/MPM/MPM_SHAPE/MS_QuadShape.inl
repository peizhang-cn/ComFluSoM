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
	// 1D quadratic B-spline shape function
	inline double ShapeQ(double x, double xc, double lx)
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
	inline double DShapeQ(double x, double xc, double lx)
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
	inline double ShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeQ(x(0),xc(0),l(0));
	}
	// 1D gradient of quadratic B-spline shape function
	inline Vector3d GradShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = DShapeQ(x(0), xc(0), l(0));
		grad(1) = 0;
		grad(2) = 0;

		return grad;
	}
	// 2D quadratic B-spline shape function
	inline double ShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1));
	}
	// 2D gradient of quadratic B-spline shape function
	inline Vector3d GradShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = ShapeQ(x(1), xc(1), l(1))*DShapeQ(x(0), xc(0), l(0));
		grad(1) = ShapeQ(x(0), xc(0), l(0))*DShapeQ(x(1), xc(1), l(1));
		grad(2) = 0;

		return grad;
	}
	// 3D quadratic B-spline shape function
	inline double ShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeQ(x(0),xc(0),l(0)) * ShapeQ(x(1),xc(1),l(1)) * ShapeQ(x(2),xc(2),l(2));
	}
	// 3D gradient of quadratic B-spline shape function
	inline Vector3d GradShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = ShapeQ(x(1), xc(1), l(1))*ShapeQ(x(2), xc(2), l(2))*DShapeQ(x(0), xc(0), l(0));
		grad(1) = ShapeQ(x(0), xc(0), l(0))*ShapeQ(x(2), xc(2), l(2))*DShapeQ(x(1), xc(1), l(1));
		grad(2) = ShapeQ(x(0), xc(0), l(0))*ShapeQ(x(1), xc(1), l(1))*DShapeQ(x(2), xc(2), l(2));

		return grad;
	}
}