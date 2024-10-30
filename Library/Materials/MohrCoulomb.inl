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
    inline void MohrCoulomb(Matrix3d& de, double Mu, double La, double Phi, double Psi, double C, Matrix3d& stress)
    {
        // Apply elastic model first
        stress += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();

        SelfAdjointEigenSolver<Matrix3d> eigensolver(stress);

        double s1 = eigensolver.eigenvalues()(2);
        double s2 = eigensolver.eigenvalues()(1);
        double s3 = eigensolver.eigenvalues()(0);

        Vector3d sb (s1, s2, s3);

        double sin0 = sin(Phi);
        double cos0 = cos(Phi);
        double sin1 = sin(Psi);
        // double cos1 = cos(Psi);

        double f = (s1-s3) +(s1+s3)*sin0 -2.*C*cos0;    

        if (f>1.e-18)
        {
            double K = La+2./3.*Mu;
            
            Matrix3d v0;
            v0.col(0) = eigensolver.eigenvectors().col(2);
            v0.col(1) = eigensolver.eigenvectors().col(1);
            v0.col(2) = eigensolver.eigenvectors().col(0);

            Vector3d sc;

            double sin01 = sin0*sin1;
            double qA0 = (8.*Mu/3.-4.*K)*sin01;
            double qA1 = Mu*(1.+sin0)*(1.+sin1);
            double qA2 = Mu*(1.-sin0)*(1.-sin1);
            double qB0 = 2.*C*cos0;

            double gsl = 0.5*(s1-s2)/(Mu*(1.+sin1));
            double gsr = 0.5*(s2-s3)/(Mu*(1.-sin1));
            double gla = 0.5*(s1+s2-2.*s3)/(Mu*(3.-sin1));
            double gra = 0.5*(2.*s1-s2-s3)/(Mu*(3.+sin1));

            double qsA = qA0-4.*Mu*(1.+sin01);
            double qsB = f;
            
            double qlA = qA0-qA1-2.*qA2;
            double qlB = 0.5*(1.+sin0)*(s1+s2)-(1.-sin0)*s3-qB0;

            double qrA = qA0-2.*qA1-qA2;
            double qrB = (1.+sin0)*s1 - 0.5*(1.-sin0)*(s2+s3) - qB0;

            double qaA = -4.*K*sin01;
            double qaB = 2.*(s1+s2+s3)/3.*sin0 - qB0;

            double minslsr = min(gsl,gsr);
            double maxlara = max(gla,gra);

            if (minslsr>0. && qsA*minslsr+qsB<0.)
            {
                double dl = -qsB/qsA;
                double ds0 = -dl*(2.*K-4.*Mu/3.)*sin1;
                sc(0) = s1+ds0-dl*(2.*Mu*(1.+sin1));
                sc(1) = s2+ds0;
                sc(2) = s3+ds0+dl*(2.*Mu*(1.-sin1));
            }
            else if (gsl>0. && gla>=gsl && qlA*gsl+qlB>=0. && qlA*gla+qlB<=0.)
            {
                // return left edge
                double dl = -qlB/qlA;
                double ds0 = dl*(4.*Mu/3.-2.*K)*sin1;
                sc(0) = sc(1) = 0.5*(s1+s2)+ds0 - dl*Mu*(1.+sin1);
                sc(2) = s3+ds0 + 2.*dl*Mu*(1.-sin1);
                // cout << "left edge" << endl;
                // abort();
            }
            else if (gsr>0. && gra>=gsr && qrA*gsr+qrB>=0. && qrA*gra+qrB<=0.)
            {
                double dl = -qrB/qrA;
                double ds0 = dl*(4.*Mu/3.-2.*K)*sin1;
                sc(0) = s1 + ds0 - 2.*dl*Mu*(1.+sin1);
                sc(1) = sc(2) = 0.5*(s2+s3) + ds0 + dl*Mu*(1.-sin1);
                // cout << "right edge" << endl;
            }
            else if (maxlara>0. && qaA*maxlara+qaB>=-1.e-24)
            {
                sc(0) = sc(1) = sc(2) = C/tan(Phi);
                // cout << "apex" << endl;
            }
            else
            {
                cout << "undefined" << endl;
                cout << s1 << " " << s2 << " " << s3 << endl;
                cout << "minslsr:" << minslsr << endl;
                cout << "qsA*minslsr+qsB: " << qsA*minslsr+qsB << endl;
                cout << "===============" << endl;
                cout << "gsl: " << gsl << endl;
                cout << "gla: " << gla << endl;
                cout << "qlA*gsl+qlB: " << qlA*gsl+qlB << endl;
                cout << "qlA*gla+qlB: " << qlA*gla+qlB << endl;
                cout << "===============" << endl;
                cout << "gsr: " << gsr << endl;
                cout << "gra: " << gra << endl;
                cout << "qrA*gsr+qrB: " << qrA*gsr+qrB << endl;
                cout << "qrA*gra+qrB: " << qrA*gra+qrB << endl; 
                cout << "===============" << endl;      
                cout << "maxlara " << maxlara << endl;
                cout << "qaA*maxlara+qaB: " << qaA*maxlara+qaB << endl;
                cout << "f: " << f << endl;
                abort();
            }

            Matrix3d sp = Matrix3d::Zero();
            sp(0,0) = sc(0);
            sp(1,1) = sc(1);
            sp(2,2) = sc(2);

            stress = v0 * sp * v0.inverse();
            double fa = (sc(0)-sc(2)) +(sc(0)+sc(2))*sin0 -2.*C*cos0;
            if (abs(fa)>1.0e-12)
            {
                cout << "f before: " << f << endl;
                cout << "f after: " << fa << endl;
                cout << sc.transpose() << endl;
                abort();            
            }
        }
    }
}