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

#ifndef __OBB_H__
#define __OBB_H__

#include "PointInsideCheck.h"

namespace OBB
{
    // Q is from lab to body
    // void FindStructGridPointWithinOBB(Vector3d& origin, Vector3d& dx, Quaterniond& Q, Vector3d& x, Vector3d& lx, vector<Vector3i>& lin)
    // {
    //     vector<Vector3i> lin0(0);
    //     lin0.reserve(1000);
    //     // lin.clear();
    //     Vector3d xr = (x-origin).array()/dx.array();
    //     Vector3d lxr = lx.array()/dx.array();
    //     Quaterniond Qi = Q.inverse();

    //     // get a AABB to surround OBB
    //     Vector3d minX (1.e30, 1.e30, 1.e30);
    //     Vector3d maxX = -minX;

    //     int ind[2] = {-1, 1};
    //     for (size_t i=0; i<=1; ++i)
    //     for (size_t j=0; j<=1; ++j)
    //     for (size_t k=0; k<=1; ++k)
    //     {
    //         Vector3d xnr = lxr;
    //         xnr(0) *= ind[i];
    //         xnr(1) *= ind[j];
    //         xnr(2) *= ind[k];

    //         Vector3d xn = Qi._transformVector(xnr);

    //         if (xn(0)<minX(0))   minX(0) = xn(0);
    //         if (xn(1)<minX(1))   minX(1) = xn(1);
    //         if (xn(2)<minX(2))   minX(2) = xn(2);

    //         if (xn(0)>maxX(0))   maxX(0) = xn(0);
    //         if (xn(1)>maxX(1))   maxX(1) = xn(1);
    //         if (xn(2)>maxX(2))   maxX(2) = xn(2);
    //     }

    //     minX += xr;
    //     maxX += xr;

    //     size_t count = 0;
    //     size_t count1 = 0;

    //     // check within
    //     for (int i=ceil(minX(0)); i<=floor(maxX(0)); ++i)
    //     for (int j=ceil(minX(1)); j<=floor(maxX(1)); ++j)
    //     for (int k=ceil(minX(2)); k<=floor(maxX(2)); ++k)
    //     {
    //         count++;
    //         Vector3d xg (i,j,k);
    //         bool inside = PointInsideCheck::PointIsInsideCuboid(xg, xr, lxr, Q);
    //         if (inside)
    //         {
    //             // count1++;
    //             // vector<int> index {i,j,k};
    //             Vector3i index (i,j,k);
    //             // lin0.emplace_back(index);
    //         }
    //     }
    //     // lin = lin0;
    // }

    void FindStructGridRangeForOBB(Vector3d& origin, Vector3d& dx, Quaterniond& Q, Vector3d& x, Vector3d& lx, Vector3i& maxX0, Vector3i& minX0)
    {
        // vector<Vector3i> lin0(0);
        // lin0.reserve(1000);
        // lin.clear();
        Vector3d xr = (x-origin).array()/dx.array();
        Vector3d lxr = lx.array()/dx.array();
        Quaterniond Qi = Q.inverse();

        // get a AABB to surround OBB
        Vector3d minX (1.e30, 1.e30, 1.e30);
        Vector3d maxX = -minX;

        int ind[2] = {-1, 1};
        for (size_t i=0; i<=1; ++i)
        for (size_t j=0; j<=1; ++j)
        for (size_t k=0; k<=1; ++k)
        {
            Vector3d xnr = lxr;
            xnr(0) *= ind[i];
            xnr(1) *= ind[j];
            xnr(2) *= ind[k];

            Vector3d xn = Qi._transformVector(xnr);

            if (xn(0)<minX(0))   minX(0) = xn(0);
            if (xn(1)<minX(1))   minX(1) = xn(1);
            if (xn(2)<minX(2))   minX(2) = xn(2);

            if (xn(0)>maxX(0))   maxX(0) = xn(0);
            if (xn(1)>maxX(1))   maxX(1) = xn(1);
            if (xn(2)>maxX(2))   maxX(2) = xn(2);
        }

        minX += xr;
        maxX += xr;

        minX0(0) = ceil(minX(0));
        minX0(1) = ceil(minX(1));
        minX0(2) = ceil(minX(2));

        maxX0(0) = ceil(maxX(0));
        maxX0(1) = ceil(maxX(1));
        maxX0(2) = ceil(maxX(2));
    }
}

#endif