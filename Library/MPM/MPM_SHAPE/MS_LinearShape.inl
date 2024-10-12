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
	// 1D linear shape function
	inline double ShapeL(double x, double xc, double lx)
	{
		double nx = 0.;
		double n = abs(x-xc)/lx;
		if (n<1.)	nx = 1.-n;
		return nx;
	}
	// 1D derivative of linear shape function
	// Not sure the derivative at zero point.
	inline double DShapeL(double x, double xc, double lx)
	{
		double d = 0;
		double n = x-xc;
		if (abs(n)<=lx)	d = -Sign(n)/lx;
		return d;
	}

	inline void LS1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
	{
		double nx = ShapeL(x(0), xc(0), l(0));
		double dx = DShapeL(x(0), xc(0), l(0));
		n = nx;
		gn(0) = dx;
		gn(1) = 0.;
		gn(2) = 0.;
	}

	inline void LS2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
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

	inline void LS3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
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
	inline double ShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeL(x(0),xc(0),l(0));
	}

	// 1D gradient of linear shape function
	inline Vector3d GradShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = DShapeL(x(0), xc(0), l(0));
		grad(1) = 0;
		grad(2) = 0;

		return grad;
	}
	// 2D linear shape function
	inline double ShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1));
	}

	// 2D gradient of linear shape function
	inline Vector3d GradShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = ShapeL(x(1), xc(1), l(1))*DShapeL(x(0), xc(0), l(0));
		grad(1) = ShapeL(x(0), xc(0), l(0))*DShapeL(x(1), xc(1), l(1));
		grad(2) = 0;

		return grad;
	}
	// 3d linear shape function
	inline double ShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		return ShapeL(x(0),xc(0),l(0)) * ShapeL(x(1),xc(1),l(1)) * ShapeL(x(2),xc(2),l(2));
	}

	// 3D gradient of linear shape function
	inline Vector3d GradShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp)
	{
		Vector3d grad;

		grad(0) = ShapeL(x(1), xc(1), l(1))*ShapeL(x(2), xc(2), l(2))*DShapeL(x(0), xc(0), l(0));
		grad(1) = ShapeL(x(0), xc(0), l(0))*ShapeL(x(2), xc(2), l(2))*DShapeL(x(1), xc(1), l(1));
		grad(2) = ShapeL(x(0), xc(0), l(0))*ShapeL(x(1), xc(1), l(1))*DShapeL(x(2), xc(2), l(2));

		return grad;
	}
}