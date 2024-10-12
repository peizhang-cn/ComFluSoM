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
	// MLS
	// polynomials
	inline VectorXd PQ2D(Vector3d& x, Vector3d& xc)
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
	inline double WeightQ(double r, double h)
	{
		double w = 0.;
		if (r<0.5*h)		w = 0.75-r*r/(h*h);
		else if (r<1.5*h)	w = 0.5*pow(1.5-r/h,2);
		return w;
	}
	// 3D quadratic spline weight function
	inline double WeightQ2D(Vector3d& x, Vector3d& xc, Vector3d& l)
	{
		double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1));
		return w;
	}
	// 3D quadratic spline weight function
	inline double WeightQ3D(Vector3d& x, Vector3d& xc, Vector3d& l)
	{
		double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1))*WeightQ(abs(x(2)-xc(2)), l(2));
		return w;
	}
	// cubic spline weight function (need mofided for h)
	inline double WeightC(double r, double h)
	{
		double w = 0.;
		if (r<0.5)		w = 0.66666666666+4.*r*r*(r-1.);
		else if (r<1.)	w = 0.66666666666+4.*r*(-1.+r-r*r/3.);
		return w;
	}
	// 3D cubic spline weight function
	inline double WeightC2D(Vector3d& x, Vector3d& xc, Vector3d& l)
	{
		double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1));
		return w;
	}
	// 3D cubic spline weight function
	inline double WeightC3D(Vector3d& x, Vector3d& xc, Vector3d& l)
	{
		double w = WeightQ(abs(x(0)-xc(0)), l(0))*WeightQ(abs(x(1)-xc(1)), l(1))*WeightQ(abs(x(2)-xc(2)), l(2));
		return w;
	}

	inline VectorXd MLS(vector<Vector3d>& Xp, Vector3d& xc, size_t wtype)
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

	inline VectorXd MLS1(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype)
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

	inline VectorXd MLS2(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype)
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
}