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
    template<typename T>
    inline void VKElastic(T& td, double mu, double la, T& stress)
    {
        T strain = 0.5*(td.transpose()*td-T::Identity());
        stress = td*(2.*mu*strain + la*strain.trace()*T::Identity())*td.transpose()/td.determinant();
    }
}

namespace Material
{
    template<typename T>
    inline void YeohModel(T& td, double C10, double C20, double C30, double K, T& stress)
    {
        T B = td*td.transpose();
        // EigenSolver<T> eigensolver(B);
        SelfAdjointEigenSolver<T> eigensolver(B);    
        // Compute strain invariants using eigenvalues
        double I1 = eigensolver.eigenvalues().squaredNorm();
        double J = td.determinant();
        double rescale = pow(J, -2./3.);
        double I1star = I1*rescale;
        double dwdi1 = C10 + 2.*C20*(I1star-3.) + 3.*C30*pow(I1star-3., 2);  // dw/dI1*
        double dwdj = K*(J-1.);                                             // dw/dJ
        T identity = td;
        identity.setIdentity();
        stress = (2./J)*dwdi1*rescale*B + (dwdj - (2.*I1star)/(3.*J)*dwdi1)*identity;
    }
}

namespace Material
{
    template<typename T>
    inline void DemirayModel(T& td, double D1, double D2, double K, T& stress)
    {
        T B = td*td.transpose();
        EigenSolver<T> eigensolver(B);  
        // Compute strain invariants using eigenvalues
        double I1 = eigensolver.eigenvalues().real().squaredNorm();
        double J = td.determinant();
        double rescale = pow(J, -2./3.);
        double I1star = I1*rescale;
        double dwdi1 = D1*D2*exp(D2*(I1star-3.));    // dw/dI1*
        double dwdj = K*(J - 1./J);                  // dw/dJ
        T identity = td;
        identity.setIdentity();
        stress = (2./J)*dwdi1*rescale*B + (dwdj - (2.*I1star)/(3.*J)*dwdi1)*identity;
        // // stress.setZero();
        // cout << "stress in" << endl;
        // cout << td << endl;
        // cout << "I1star: " << I1star << endl;
        // cout << "J: " << J << endl;
        // cout << "dwdi: " << dwdi << endl;
        // cout << stress << endl;
        // cout << "===========" << endl;
    }
}