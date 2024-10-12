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

class MPM_CELL
{
public:
    MPM_CELL();
    ~MPM_CELL();
    MPM_CELL(const Vector3d& x);
    void Reset();

    size_t                      ID;                     // Index of gird in the list

    size_t                      Np;

    double                      Vol0;                   // Inital volume mapped from particles
    double                      Vol;                    // Current volume mapped from particles
    double                      J;                      // Determinant of gradient deformation

    Vector3d                    X;                      // Cell position

    Matrix3d                    Stress;
};

inline MPM_CELL::MPM_CELL(const Vector3d& x)
{
    ID      = 0;
    Np      = 0;
    Vol0    = 0.;
    Vol     = 0.;
    X       = x;
    Stress.setZero();
}

inline void ::MPM_CELL::Reset()
{
    Np   = 0;
    Vol0 = 0.;
    Vol  = 0.;
    J    = 1.;
    Stress.setZero();
}

