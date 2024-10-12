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
	// 1D GIMP shape function
	inline double ShapeGIMP(double x, double xc, double lx, double lpx)
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
	inline double DShapeGIMP(double x, double xc, double lx, double lpx)
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
	inline void GIMP(double x, double xc, double lx, double lpx, double& n, double& gn)
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

	inline void GIMP1D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
	{
		double n0 = 0.;
		double gn0 = 0.;

		GIMP(x(0), xc(0), l(0), lp(0), n0, gn0);

		n = n0;

		gn(0) = gn0;
		gn(1) = 0.;
		gn(2) = 0.;
	}

	inline void GIMP2D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
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

	inline void GIMP3D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn)
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
}