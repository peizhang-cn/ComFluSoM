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

#include <DEM.h>

int main(int argc, char const *argv[])
{
    double lx = 0.02;                           // domain size x (m)
    double ly = 0.02;                           // domain size y (m)
    double lz = 0.05;                           // domain size z (m)
    double dt = 1.e-6;                          // time step (s)
    double cr = 0.2;                            // restitution coefficient
    double miu = 0.2;                           // friction coefficient
    double r = 1.e-3;                           // radius of particle
    double rho = 2000.;                         // density (kg/m^3)

    Vector3d origin (0.,0.,0.);                 // starting point of domain
    Vector3d domSize (lx, ly, lz);              // domain size
    // use Linear contact model
    DEM* dem = new DEM(origin, domSize, "LINEAR_CR", dt);
    dem->SetMaterialCoef(0, 0, cr, miu, miu);   // set materials
    dem->SetPeriodic(false, false, false);      // remove periodic boundary
    int nx = floor(0.5*lx/r)-1;
    int ny = floor(0.5*ly/r)-1;
    int nz = floor(0.5*lz/r)-1;
    Vector3i binSize (nx,ny,nz);                // size of bin for contact detection
    dem->InitBinSystem(binSize);                // generate bin system
    cout << "bin dimension: " << dem->Bins->Dx.transpose() << endl;
    dem->SetParallel(8);                        // choose how many processors to use
    // z-min wall
    Vector3d xcube0 = 0.5*domSize + origin;
    xcube0(2) = origin(2);
    Vector3d lcube0 = domSize;
    lcube0(2) = 2.*r;
    dem->AddCuboid(-2, lcube0, xcube0, rho);
    dem->Lp[dem->Lp.size()-1]->Fix();
    // x-min wall
    Vector3d xcube1 = 0.5*domSize + origin;
    xcube1(0) = origin(0);
    Vector3d lcube1 = domSize;
    lcube1(0) = 2.*r;
    dem->AddCuboid(-3, lcube1, xcube1, rho);
    dem->Lp[dem->Lp.size()-1]->Fix();
    // x-max wall
    Vector3d xcube2 = 0.5*domSize + origin;
    xcube2(0) = origin(0)+domSize(0);
    Vector3d lcube2 = domSize;
    lcube2(0) = 2.*r;
    dem->AddCuboid(-4, lcube2, xcube2, rho);
    dem->Lp[dem->Lp.size()-1]->Fix();
    // y-min wall
    Vector3d xcube3 = 0.5*domSize + origin;
    xcube3(1) = origin(1);
    Vector3d lcube3 = domSize;
    lcube3(1) = 2.*r;
    dem->AddCuboid(-5, lcube3, xcube3, rho);
    dem->Lp[dem->Lp.size()-1]->Fix();
    // y-max wall
    Vector3d xcube4 = 0.5*domSize + origin;
    xcube4(1) = origin(1)+domSize(1);
    Vector3d lcube4 = domSize;
    lcube4(1) = 2.*r;
    dem->AddCuboid(-6, lcube4, xcube4, rho);
    dem->Lp[dem->Lp.size()-1]->Fix();

    // the generated packing is not a dense packing, so the box should be large enough
    size_t np = 500;                            // number of particle
    Vector3d x0 (r, r, r);                      // minimum position of box
    Vector3d x1 (lx-r,ly-r,0.5*lz);             // maximum position of box
    dem->AddNSpheres(-1, 0, np, x0, x1, r, 0., rho);
    // set gravity
    double g = 9.8;
    Vector3d gravity (0.,0.,-g);
    dem->SetG(gravity);
    // set contact parameters
    for (size_t p=0; p<dem->Lp.size(); ++p)
    {
        dem->Lp[p]->Kn = 1.e6;                  // normal stiffness
        dem->Lp[p]->Kt = 1.e2;                  // tangential stiffness
    }

    // solve ========================================================================
    size_t tt = 4e5;                            // total time steps
    size_t ts = 2e3;                            // saved time steps
    for (size_t t=0; t<tt; ++t)
    {
        if (t%ts==0)    cout << "time step: " << t << endl;
        vector<double> times(0);
        dem->SolveOneStep(t, ts, times);
    }

  return 0;
}
