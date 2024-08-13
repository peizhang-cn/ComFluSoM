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

// this is a test program for DEM
// it shows a sphere moving on a plane with initial velocity and friction
// the simulation consist two phase: 
// 1 the sphere is falling down under gravity
// 2 the sphere is sliding on the plane with zero angular velocity
// at the end of the simulation, the thoerictical displacement is compared with the simulation result

#include <DEM.h>

int main(int argc, char const *argv[])
{
    double lx = 0.02;                           // domain size x (m)
    double ly = 0.002;                          // domain size y (m)
    double lz = 0.02;                           // domain size z (m)
    double dt = 1.e-6;                          // time step (s)
    double cr = 0.2;                            // restitution coefficient
    double miu = 0.2;                           // friction coefficient

    Vector3d origin (0.,0.,0.);                 // starting point of domain
    Vector3d domSize (lx, ly, lz);              // domain size
    // use Linear contact model
    DEM* dem = new DEM(origin, domSize, "LINEAR_CR", dt);
    dem->SetMaterialCoef(0, 0, cr, miu, miu);   // set materials
    dem->SetPeriodic(false, false, false);      // remove periodic boundary
    Vector3i binSize (1,1,1);                   // size of bin for contact detection
    dem->InitBinSystem(binSize);                // generate bin system
    dem->SetParallel(1);                        // choose how many processors to use
    // add a sphere
    double rho = 2000.;                         // density (kg/m^3)
    double r = 1.e-3;                           // radius of particle
    Vector3d xp0 (0.05*lx, 0.5*ly, 2.*r);       // mass center of particle
    dem->AddSphere(-1/* tag */, r, xp0, rho);   // add a sphere 
    // add a cuboid
    Vector3d xcube (0.5*lx, 0.5*ly, 0.);        // mass center of cube
    double width = 0.2*r;                       // width of cube
    Vector3d lcube (lx, ly, width);             // length of cube
    dem->AddCuboid(-2, lcube, xcube, rho);      // add a cuboid
    dem->Lp[dem->Lp.size()-1]->Fix();           // fix the cuboid
    // set gravity
    double g = 9.8;
    Vector3d gravity (0.,0.,-g);
    dem->SetG(gravity);
    // set contact parameters
    for (size_t p=0; p<dem->Lp.size(); ++p)
    {
        dem->Lp[p]->Kn = 1.e7;                  // normal stiffness
        dem->Lp[p]->Kt = 1.e3;                  // tangential stiffness
    }

    // solve ========================================================================
    size_t tt = 1e5;                            // total time steps
    size_t ts = 2e5;                            // saved time steps
    // main loop
    for (size_t t=0; t<tt; ++t)
    {
        if (t%ts==0)    cout << "first solve time step: " << t << endl;
        vector<double> times(0);
        dem->SolveOneStep(t, ts, times);
    }
    // reset material properties to remove dampping
    cr = 1.;
    dem->SetMaterialCoef(0, 0, cr, miu, miu);   // set materials
    double v0 = 0.15;                           // initial velocity
    dem->Lp[0]->V(0) = v0;                      // set initial velocity
    double x0 = dem->Lp[0]->X(0);               // initial position
    Vector3d w0 (0.,0.,0.);
    dem->Lp[0]->FixW(w0);                       // fix particle's rotation
    // reset time steps
    tt = 1e5;                                   // total time steps
    ts = 2000;                                  // saved time steps
    for (size_t t=0; t<tt; ++t)
    {
        if (t%ts==0)    cout << "second solve time step: " << t << endl;
        vector<double> times(0);
        dem->SolveOneStep(t, ts, times);
    }
    double xf = dem->Lp[0]->X(0);               // final position
    double dis_t = 0.5*v0*v0/(miu*g);           // theoretical displacement
    double dis_s = xf-x0;                       // simulated displacement
    cout << "\033[31m" << "thoerictical displacement: " << "\033[0m" << dis_t << "\033[31m" << endl;
    cout << "\033[31m" << "simulated displacement: " << "\033[0m" << dis_s << "\033[31m" << endl;
    cout << "\033[31m" << "error:" << "\033[0m" << abs(dis_s/dis_t-1.) << endl;
    return 0;
}
