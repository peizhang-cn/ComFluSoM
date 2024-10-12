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
 
#pragma once

namespace MPM_ShapeFunction
{
	double Sign(double x)
	{
		return (x<0.)? -1.:1.;
	}

	double ShapeL(double x, double xc, double lx);
	double DShapeL(double x, double xc, double lx);
	void LS1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);
	void LS2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);
	void LS3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);
	double ShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeL1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	double ShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeL2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	double ShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeL3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);

	double ShapeQ(double x, double xc, double lx);
	double DShapeQ(double x, double xc, double lx);
	double ShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeQ1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	double ShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeQ2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	double ShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d GradShapeQ3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);

	double ShapeCubicBSpline1(double x, double xc, double lx);
	double ShapeCubicBSpline2(double x, double xc, double lx);
	double ShapeCubicBSpline3(double x, double xc, double lx);
	double ShapeCubicBSpline4(double x, double xc, double lx);
	double ShapeCubicBSpline(int type, double x, double xc, double lx);

	double DShapeCubicBSpline1(double x, double xc, double lx);
	double DShapeCubicBSpline2(double x, double xc, double lx);
	double DShapeCubicBSpline3(double x, double xc, double lx);
	double DShapeCubicBSpline4(double x, double xc, double lx);
	double DShapeCubicBSpline(int type, double x, double xc, double lx);

	double ShapeCubicBSpline1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	Vector3d GradShapeCubicBSpline1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	void CubicBSpline1D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn);

	double ShapeCubicBSpline2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	Vector3d GradShapeCubicBSpline2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	void CubicBSpline2D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn);

	double ShapeCubicBSpline3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	Vector3d GradShapeCubicBSpline3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type);
	void CubicBSpline3D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& type, double& n, Vector3d& gn);	
	// double ShapeC(double x, double xc, double lx);
	// double DShapeC(double x, double xc, double lx);
	// double ShapeC1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	// Vector3d GradShapeC1D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	// double ShapeC2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	// Vector3d GradShapeC2D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	// double ShapeC3D (Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);

	double ShapeGIMP(double x, double xc, double lx, double lpx);
	double DShapeGIMP(double x, double xc, double lx, double lpx);
	void GIMP(double x, double xc, double lx, double lpx, double& n, double& gn);
	void GIMP1D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);
	void GIMP2D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);
	void GIMP3D(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);

	VectorXd PQ2D(Vector3d& x, Vector3d& xc);
	double WeightQ(double r, double h);
	double WeightQ2D(Vector3d& x, Vector3d& xc, Vector3d& l);
	double WeightQ3D(Vector3d& x, Vector3d& xc, Vector3d& l);
	double WeightC(double r, double h);
	double WeightC2D(Vector3d& x, Vector3d& xc, Vector3d& l);
	double WeightC3D(Vector3d& x, Vector3d& xc, Vector3d& l);
	VectorXd MLS(vector<Vector3d>& Xp, Vector3d& xc, size_t wtype);
	VectorXd MLS1(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype);
	VectorXd MLS2(vector<Vector3d>& Xp, Vector3d& xc, double vis, size_t wtype);
}

#include "MS_LinearShape.inl"
#include "MS_QuadShape.inl"
#include "MS_CubicShape.inl"
#include "MS_GimpShape.inl"
#include "MS_MlsShape.inl"