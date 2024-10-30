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

namespace Material
{
    inline void DruckerPrager(Matrix3d& de, double Mu, double La, double A_dp, double B_dp, double Ad_dp, double C, Matrix3d& stress)
    {
        // Apply elastic model first
        stress += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();

        double p = stress.trace()/3.;
        Matrix3d ss = stress - p*Matrix3d::Identity();
        double j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
        double f = j2sqr+A_dp*p-B_dp*C;
        if (f>0.)
        {
            double K = La+2./3.*Mu;
            // return to cone
            if (A_dp*(p-j2sqr/Mu*K*Ad_dp)-B_dp*C<0.)
            {
                double dl = (j2sqr+A_dp*p-B_dp*C)/(Mu+A_dp*K*Ad_dp);
                stress -= dl*(Mu/j2sqr*ss+K*Ad_dp*Matrix3d::Identity());
            }
            // return to apex
            else
            {
                stress = B_dp*C/A_dp*Matrix3d::Identity();
            }
            p = stress.trace()/3.;
            ss = stress - p*Matrix3d::Identity();
            j2sqr = sqrt(0.5*(ss.array()*ss.array()).sum());
            double fb = f;
            f = j2sqr+A_dp*p-B_dp*C;
            if (f>1.0e-8)
            {
                cout << "f before: " << fb << endl;
                cout << "f after: " << f << endl;
                abort();
            }
        }
    }
}