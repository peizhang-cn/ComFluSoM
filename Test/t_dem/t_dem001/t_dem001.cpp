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
// it shows two spheres contact with opposite initial velocity
// the engergy loss is controlled by the coefficient of restitution

#include <DEM.h>

int main(int argc, char const *argv[])
{
    double lx = 0.02;                           // domain size x (m)
    double ly = 0.02;                           // domain size y (m)
    double lz = 0.02;                           // domain size z (m)
    double dt = 1.e-7;                          // time step (s)
    double cr = 0.2;                            // restitution coefficient
    double miu = 0.0;                           // friction coefficient

    Vector3d origin (0.,0.,0.);                 // starting point of domain
    Vector3d domSize (lx, ly, lz);              // domain size
    // use Linear contact model
    DEM* dem = new DEM(origin, domSize, "LINEAR_CR", dt);
    dem->SetMaterialCoef(0, 0, cr, miu, miu);   // set materials
    dem->SetPeriodic(false, false, false);      // remove periodic boundary
    Vector3i binSize (1,1,1);                   // size of bin for contact detection
    dem->InitBinSystem(binSize);                // generate bin system
    dem->SetParallel(1);                        // choose how many processors to use
    // add first sphere
    double rho = 2000.;                         // density (kg/m^3)
    double r = 1.e-3;                           // radius of particle
    Vector3d xp0 (0.5*lx-1.1*r, 0.5*ly, 0.5*lz);// mass center of particle
    dem->AddSphere(-1/* tag */, r, xp0, rho);   // add a sphere 
    // add second sphere
    Vector3d xp1 (0.5*lx+1.1*r, 0.5*ly, 0.5*lz);// mass center of particle
    dem->AddSphere(-2/* tag */, r, xp1, rho);   // add a sphere 
    double v0 = 0.1;                            // init velocity
    dem->Lp[0]->V(0) = v0;                      // set velocity
    dem->Lp[1]->V(0) = -v0;                     // set velocity
    // set contact parameters
    for (size_t p=0; p<dem->Lp.size(); ++p)
    {
        dem->Lp[p]->Kn = 1.e7;                  // normal stiffness
        dem->Lp[p]->Kt = 1.e3;                  // tangential stiffness
    }

    // solve ========================================================================
    size_t tt = 1e5;                            // total time steps
    size_t ts = 5000;                           // saved time steps
    // main loop
    for (size_t t=0; t<tt; ++t)
    {
        if (t%ts==0)    cout << "time step: " << t << endl;
        vector<double> times(0);
        dem->SolveOneStep(t, ts, times);
    }
    double cr_sim = abs(dem->Lp[0]->V(0))/v0;   // simulated restitution coefficient
    cout << "\033[31m" << "restitution coefficient: " << "\033[0m" << cr << endl;
    cout << "\033[31m" << "simulated restitution coefficient: " << "\033[0m" << cr_sim << endl;
    cout << "\033[31m" << "error:" << "\033[0m" << abs(cr_sim/cr-1.) << endl;
    return 0;
}
