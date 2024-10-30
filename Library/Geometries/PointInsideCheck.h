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

#include "Intersection.h"

namespace PointInsideCheck
{
    template <typename T>
    bool Sphere(T x, T xc, double r)
    {
        bool inside = true;
        if ((x-xc).norm()>r)     inside = true;
        return inside;
    }

    // under body frame, h is half length
    template <typename T>
    bool Cylinder(T xr, double h, double r)
    {
        bool inside = false;
        double hr = xr(2);
        double r2r = xr(0)*xr(0) + xr(1)*xr(1);
        if (abs(hr)<=h && r2r<=r*r) inside = true;
        return inside;
    }

    // under body frame, l is half length
    template <typename T>
    bool Cuboid(T xr, T l)
    {
        bool inside = true;
        for (int d=0; d<l.size(); ++d)
        {
            double dis = l(d)-abs(xr(d));
            if (dis<0.)
            {
                inside = false;
                break;
            }
        }
        return inside;
    }

    template <typename T>
    Vector3d BarycentricCoordinates2D(T x, T p0, T p1, T p2)
    {
        Vector2d xm (x(0),x(1));
        Vector2d pm0 (p0(0),p0(1));
        Vector2d pm1 (p1(0),p1(1));
        Vector2d pm2 (p2(0),p2(1));

        Matrix2d A;
        A.col(0) = pm0-pm2;
        A.col(1) = pm1-pm2;
        Vector2d b = xm-pm2;
        // Vector2d u2 = A.colPivHouseholderQr().solve(b); not accurate
        Vector2d u2 = A.inverse()*b;
        Vector3d u (u2(0),u2(1), 1.-u2(0)-u2(1));
        return u;
    }

    template <typename T>
    Vector3d Barycentric(T p, T a, T b, T c)
    {
        T v0 = b - a;
        T v1 = c - a;
        T v2 = p - a;
        double d00 = v0.dot(v0);
        double d01 = v0.dot(v1);
        double d11 = v1.dot(v1);
        double d20 = v2.dot(v0);
        double d21 = v2.dot(v1);
        double invDenom  = 1./(d00 * d11 - d01 * d01);

        Vector3d u;
        u(1) = (d11 * d20 - d01 * d21) * invDenom;
        u(2) = (d00 * d21 - d01 * d20) * invDenom;
        u(0) = 1. - u(1) - u(2);
        return u;
    }

    template <typename T>
    bool Triangle(T x, T p0, T p1, T p2)
    {
        bool isInside = true;
        Vector3d u = BarycentricCoordinates2D(x, p0, p1, p2);
        if (u(0)<0. || u(1)<0. || u(2)<0.)      isInside = false;
        return isInside;
    }

    // check if a point is inside a 2D quadrilateral
    template <typename T>
    bool Quadrilateral(T x, T p0, T p1, T p2, T p3)
    {
        bool isInside0 = Triangle(x, p0, p1, p2);
        bool isInside1 = Triangle(x, p2, p3, p0);
        return isInside0&&isInside1;
    }

    // l is half length of the cuboid
    template <typename T>
    bool Cuboid(T x, T xc, T l, Quaterniond quat)
    {
        bool isInside = true;
        T xr = quat._transformVector(x-xc);
        for (int d=0; d<l.size(); ++d)
        {
            if (abs(xr(d))>l(d))
            {
                isInside = false;
                break;
            }
        }
        return isInside;
    }

    // check if a point is inside a tetrahedron
    // https://math.stackexchange.com/questions/3698021/how-to-find-if-a-3d-point-is-in-on-outside-of-tetrahedron
    template <typename T>
    bool Tetrahedron(T x, T p0, T p1, T p2, T p3)
    {
        bool isInside = true;
        Matrix3d A;
        A.col(0) = p1-p0;
        A.col(1) = p2-p0;
        A.col(2) = p3-p0;
        T b = x-p0;
        // Vector3d u = A.colPivHouseholderQr().solve(b);
        T u = A.inverse()*b;
        if (u(0)<0. || u(1)<0. || u(2)<0. || u(0)+u(1)+u(2)>1.)     isInside = false;
        return isInside;
    }

    template <typename T>
    bool ConvexPolygon2D(T& x, vector<T>& P, vector<Vector2i>& E)
    {
        bool isInside = true;
        for (size_t e=0; e<E.size(); ++e)
        {
            T v0 = P[E[e](1)] - P[E[e](0)];
            T v1 = x - P[E[e](0)];
            double sign = v1(0)*v0(1) - v1(1)*v0(0);
            if (sign>0.)
            {
                isInside = false;
                break;
            }
        }
        return isInside;
    }

    // need to closer look to avoid corner case
    // template <typename T>
    // bool PointIsInsidePolygon2D(T x0, double ang, vector<T> P, vector<Vector2i> E)
    // {
    //     bool isInside = false;
    //     size_t crossNumLeft = 0;
    //     size_t crossNumRight = 0;

    //     Vector2d x (x0(0),x0(1));

    //     Rotation2Dd rot(ang);

    //     for (size_t e=0; e<E.size(); ++e)
    //     {
    //         Vector2d c0 (P[E[e](0)](0)+1.e-14, P[E[e](0)](1)+1.e-14);
    //         Vector2d d0 (P[E[e](1)](0)+1.e-14, P[E[e](1)](1)+1.e-14);

    //         Vector2d c = rot.toRotationMatrix()*(c0-x)+x;
    //         Vector2d d = rot.toRotationMatrix()*(d0-x)+x;

    //         double s = (x(1)-d(1))/(c(1)-d(1));
    //         double t = s*(c(0)-d(0))+d(0);
    //         if (s>0. && s<=1.)
    //         {
    //             if (t<x(0)) crossNumLeft++;
    //             else        crossNumRight++;
    //         }
    //     }
    //     if (crossNumLeft%2==1 && crossNumRight%2==1)    isInside = true;
    //     return isInside;
    // }

    // not tested yet!!!
    template <typename T>
    bool Polygon2D(T x0, double ang, vector<T> P, vector<Vector2i> E)
    {
        bool isInside = false;
        size_t crossNumLeft = 0;
        size_t crossNumRight = 0;

        Vector2d x (x0(0),x0(1));

        Vector2d xe (x(0)+cos(ang), x(1)+sin(ang));

        for (size_t e=0; e<E.size(); ++e)
        {
            Vector2d c0 (P[E[e](0)](0), P[E[e](0)](1));
            Vector2d d0 (P[E[e](1)](0), P[E[e](1)](1));

            Vector2d u = Intersection::LineLine2D(x,xe,c0,d0);

            if (u(1)>0. && u(1)<=1.)
            {
                if (u(1)<0.) crossNumLeft++;
                else        crossNumRight++;
            }
        }
        if (crossNumLeft%2==1 && crossNumRight%2==1)    isInside = true;
        return isInside;
    }

    // need to improve efficiency
    template <typename T>
    bool Polygon2D(T x, vector<T> P, vector<Vector2i> E)
    {
        // try five times then vote
        bool isInside = false;
        bool isInside0 = Polygon2D(x,0.00*M_PI,P,E);
        bool isInside1 = Polygon2D(x,0.5*M_PI,P,E);
        bool isInside2 = Polygon2D(x,0.75*M_PI,P,E);
        bool isInside3 = Polygon2D(x,1.25*M_PI,P,E);
        bool isInside4 = Polygon2D(x,1.75*M_PI,P,E);
        // bool isInside5 = PointIsInsidePolygon2D(x,0.5 *M_PI,P,E);

        size_t num = isInside0+isInside1+isInside2+isInside3+isInside4;
        if (num>2)  isInside = true;

        return isInside;
    }

    // Warning: may not work while the point is on surface
    template <typename T>
    bool Polyhedron(T x, vector<T> P, vector<VectorXi> F)
    {
        bool isInside = false;
        double sum = 0.;
        for (size_t f=0; f<F.size(); ++f)
        {
            Vector3d A = P[F[f](0)]-x;
            Vector3d B = P[F[f](1)]-x;
            Vector3d C = P[F[f](2)]-x;

            double a = A.norm();
            double b = B.norm();
            double c = C.norm();

            Matrix3d mat;
            mat.col(0) = A;
            mat.col(1) = B;
            mat.col(2) = C;

            double y = mat.determinant();
            double x = a*b*c + A.dot(B)*c + B.dot(C)*a + C.dot(A)*b;
            sum += atan2(y,x);
        }
        if (abs(sum)>=1.e-3)    isInside = true;
        return isInside;
    }
}