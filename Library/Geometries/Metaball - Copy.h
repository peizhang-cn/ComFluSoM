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

namespace MetaballFunctions
{
    double CalMetaC(vector<Vector3d>& P, VectorXd& K, Vector3d& X)
    {
        Vector3d W (1.,1.,1.);
        double c = 0.;
        for (size_t i=0; i<P.size(); ++i)
        {
            Vector3d xpi = X-P[i];
            Vector3d sqr = xpi.array()*xpi.array();
            double ri2 = W.dot(sqr);
            c += K(i)/ri2;
        }
        return c;
    }
    void CalMetaCF(vector<Vector3d>& P, VectorXd& K, Vector3d& X, double& C, Vector3d& F)
    {
        double c = 0.;
        Vector3d f = Vector3d::Zero();
        for (size_t i=0; i<P.size(); ++i)
        {
            Vector3d a = X-P[i];
            double ri2inv = 1./a.squaredNorm();
            double ri4inv = ri2inv*ri2inv;
            f -= 2.*K(i)*ri4inv*a;
            c += K(i)*ri2inv;
        }
        C = c;
        F = f;
    }

    void FindMetaVectorCrossPointsWithRay(vector<Vector3d>& P, VectorXd& K, Vector3d& x0, Vector3d& vec, Vector3d& cp, size_t inum, double tol)
    {
        // x0: starting point of the vector, vec: direction, sp: surface point
        Vector3d xi = x0;
        Vector3d xj = x0+vec;

        double c0 = CalMetaC(P, K, xi);
        double c1 = CalMetaC(P, K, xj);
        if (c0<1.)
        {
            cout << "x0 is out of the metaball!" << endl;
            abort();
        }
        // move xi and xj to make sure that ci>1 and cj<1.
        for (size_t s=1; s<1e5; ++s)
        {
            if (c1<1.)	break;
            else
            {
                xi = xj;
                xj = x0 + pow(2,s)*vec;
                c1 = CalMetaC(P, K, xj);
            }
        }
        // use middle point method to find crossing point
        for (size_t s=0; s<inum; ++s)
        {
            Vector3d xm = 0.5*(xi+xj);
            double cm = CalMetaC(P, K, xm);
            if (abs(cm-1.)<tol)
            {
                cp = xm;
                break;
            }
            if (cm>1.)	xi = xm;
            else		xj = xm;
        }
    }

    Vector3d ProjectPointToMetaball(vector<Vector3d>& P, VectorXd& K, Vector3d& x0, size_t inum, double tol)
    {
        Vector3d xp = x0;
        for (size_t s=0; s<inum; ++s)
        {
            double c; Vector3d f;
            CalMetaCF(P, K, xp, c, f);
            if (abs(c-1.)<tol)  break;
            xp += (1.-c)/f.squaredNorm()*f;
        }
        return xp;
    }

    void CalSurfaceMesh2D(vector<Vector3d>& P, VectorXd& K, double dis, vector<Vector3d>& ver)
    {
        Vector3d xc = Vector3d::Zero();
        for (size_t i=0; i<P.size(); ++i)   xc += P[i];
        xc /= P.size();
        double rMax = -1.;
        double kMax = -1.;
        for (size_t i=0; i<P.size(); ++i)
        {
            double r = (P[i]-xc).norm();
            if (r>rMax)	rMax = r;
            if (K(i)>kMax)	kMax = K(i);
        }
        double lx = 2.*(rMax + sqrt(kMax));
        double cinit = 0.;
        Vector3d pinit;
        srand( (unsigned)time(NULL) );
        while (abs(cinit-1.)>0.2)
        {
            pinit(0) = lx*(rand()/double(RAND_MAX) - 0.5);
            pinit(1) = lx*(rand()/double(RAND_MAX) - 0.5);
            pinit(2) = 0.;
            pinit += xc;
            cinit = CalMetaC(P, K, pinit);
        }


        // cout << "pinit: " << pinit.transpose() << endl;
        Vector3d sp0 = ProjectPointToMetaball(P, K, pinit, 1000, 1e-3);
        // cout << "sp0: " << sp0.transpose() << endl;
        ver.push_back(sp0);
        auto findNewPoint = [](vector<Vector3d>& P, VectorXd& K, double dis, Vector3d& x0) -> Vector3d
        {
            Vector3d zAxis (0.,0.,1.);
            double c; Vector3d f;
            CalMetaCF(P, K, x0, c, f);
            Vector3d ns = f.cross(zAxis);
            ns.normalize();
            Vector3d x1 = x0+ns*dis;
            Vector3d xp = ProjectPointToMetaball(P, K, x1, 1000, 1e-3);
            return xp;
        };
        Vector3d sp1 = findNewPoint(P, K, dis, sp0);
        ver.push_back(sp1);
        bool endFlag = false;
        while (!endFlag)
        {
            Vector3d spi = findNewPoint(P, K, dis, ver[ver.size()-1]);
            double dis0 = (spi-sp0).norm();
            double dis1 = (spi-sp1).norm();
            if (dis0<0.5*dis || dis1<0.5*dis)
            {
                endFlag = true;
                // cout << "spi: " << spi.transpose() << endl;
                // cout << "sp0: " << sp0.transpose() << endl;
                // cout << "dis0: " << dis0 << ", dis1: " << dis1 << endl;
                // abort();
            }
            else                                ver.push_back(spi);
        }
        // reverse(ver.begin(), ver.end());
    }

    void FindClosestPointFromSurfacePoints(vector<Vector3d>& P0, VectorXd& K0, vector<Vector3d>& Sp, Vector3d& Cp)
    {
        size_t id = 0;
        double cmax = 0.;
        for (size_t i=0; i<Sp.size(); ++i)
        {
            double c0 = CalMetaC(P0, K0, Sp[i]);
            if (c0>cmax)
            {
                cmax = c0;
                id = i;
            }
        }
        Cp = Sp[id];
    }

    void FindClosestPointFromPointToMetaball(Vector3d& p0, vector<Vector3d>& metaPi, VectorXd& metaKi, double step0, Vector3d& initP, Vector3d& cp, bool& convFlag)
    {
        convFlag = false;
        Vector3d cpi = initP;

        size_t nk = 2000;
        double step = step0;
        double resb = 2.;
        Vector3d cpib = cpi;
        for (size_t k=0; k<nk; ++k)
        {
            double ci; Vector3d fi;
            MetaballFunctions::CalMetaCF(metaPi, metaKi, cpi, ci, fi);

            double dxn = (1.-ci)/fi.norm();
            if (fi.norm()<1.e-12)   dxn = 0.;
            fi.normalize();           
            Vector3d fj = p0-cpi;
            fj.normalize();
            auto fjt = fj - fj.dot(fi)*fi;
            double res = fjt.norm();
                // cout << "k: " << k << endl;
                // cout << "fi: " << fi.transpose() << endl;
                // cout << "fj: " << fj.transpose() << endl;  
                // cout << "fjt.norm(): " << fjt.norm() << endl;
                // cout << "abs(ci-1.): " << abs(ci-1.) << endl;
                // cout << "=====================" << endl; 

            if (res<1e-8 && abs(ci-1.)<1.e-8)
            {
                // cout << "initP: " << initP.transpose() << endl;
                // cout << "fi: " << fi.transpose() << endl;
                // cout << "fj: " << fj.transpose() << endl;   
                // cout << "k: " << k << endl;
                convFlag = true;
                break;
            }

            cpi += step*fjt;
            cpi += dxn*fi;

            // if (res>resb)
            // {
            //     // cout << "res>resb" << endl;
            //     cpi = cpib;
            //     res = resb;
            //     step *= 0.5;
            // }
            resb = res;
            cpib = cpi;
        }
        // if (!convFlag)
        // {
        //     cout << "FindClosestPointFromPointToMetaball not converge!!!" << endl;
        // }
        cp = cpi;     
    }

    // void FindClosestPointFromPointToMetaball(Vector3d& p0, vector<Vector3d>& metaPi, VectorXd& metaKi, double step0, Vector3d& initP, Vector3d& cp, bool& convFlag)
    // {
    //     ProjectPointToMetaball(metaPi, metaKi, p0, 200, 1.e-6);
    //     convFlag = tr
    // }

    void FindClosestPointFromMetaballToMetaball(vector<Vector3d>& metaPi, vector<Vector3d>& metaPj, VectorXd& metaKi, VectorXd& metaKj, double step, Vector3d& initP, Vector3d& cp, bool& convFlag)
    {
        convFlag = false;
        Vector3d cpi = initP;

        size_t nk = 200;
        for (size_t k=0; k<nk; ++k)
        {
            double ci; Vector3d fi;
            MetaballFunctions::CalMetaCF(metaPi, metaKi, cpi, ci, fi);
            double dxn = (1.-ci)/fi.norm();
            if (fi.norm()<1.e-12)   dxn = 0.;
            fi.normalize();
            double cj; Vector3d fj;
            MetaballFunctions::CalMetaCF(metaPj, metaKj, cpi, cj, fj);
            fj.normalize();
            auto fjt = fj - fj.dot(fi)*fi;

            if (fjt.norm()<1e-5 && abs(ci-1.)<1.e-5)
            {
                // cout << "k: " << k << endl;
                convFlag = true;
                break;
            }

            cpi += step*fjt;
            cpi += dxn*fi;
        }
        cp = cpi;
    }

    // p0 must be counter-clock order, polygon must be convex, no intersections
    void FindClosestPointFromPolygonToMetaball2D(vector<Vector3d>& p0, vector<Vector3d>& metaP, VectorXd& metaK, double stepRatio, Vector3d& initP, Vector3d& cp, bool& convFlag)
    {
        double cb; Vector3d fb;
        CalMetaCF(metaP, metaK, p0[0], cb, fb);
        size_t maxID = 0;
        double maxC = cb;
        bool isLine = false;
        for (size_t i=0; i<p0.size(); ++i)
        {
            size_t id = (i+1)%p0.size();
            double ci; Vector3d fi;
            CalMetaCF(metaP, metaK, p0[id], ci, fi);
            Vector3d axisz (0.,0.,1.);
            Vector3d vec = p0[id]-p0[i];
            Vector3d norml = vec.cross(axisz);
            size_t idn = (i+2)%p0.size();
            double sideCheck =(p0[idn]-p0[i]).dot(norml)*(metaP[0]-p0[i]).dot(norml);
            if (sideCheck<0.)
            {
                // cout << "sideCheck pass" << endl;
                // cout << "i: " << i << endl;
                // cout << "norml: " << norml.transpose() << endl;
                // cout << "fb: " << fb.transpose() << endl;
                // cout << "fi: " << fi.transpose() << endl;
                // cout << "(norml.cross(fb))(2): " << (norml.cross(fb))(2) << endl;
                // cout << "(norml.cross(fi))(2): " << (norml.cross(fi))(2) << endl;
                // if ((norml.cross(fb))(2)*(norml.cross(fi))(2)<0.)

                if ((norml.cross(fb))(2)*(norml.cross(fi))(2)<0.)
                // if (cm>cb && cm>ci)
                {
                    isLine = true;
                    Vector3d vecn = vec.normalized();
                    // cout << "isline" << endl;
                    // Vector3d cp0 = (cb*p0[i] + ci*p0[id])/(ci+cb);
                    Vector3d cp0 = ClosestPoint::SegmentClosestPoint(initP, p0[i], p0[id]);
                    convFlag = false;
                    size_t convk = 0;
                    double res = 1.e20;
                    for (size_t k=0; k<2000; ++k)
                    {
                        double ck; Vector3d fk;
                        CalMetaCF(metaP, metaK, cp0, ck, fk);
                        fk.normalize();
                        res = fk.dot(vecn);
                        cp0 += res*stepRatio*vec;
                        // cout << "k: " << k << endl;
                        // cout << "cp0: " << cp0.transpose() << endl;
                        // cout << "ck: " << ck << endl;
                        // cout << "res: " << res << endl;
                        // cout << "=======" << endl;
                        if (abs(res)<1.e-6)
                        {
                            convk = k;
                            // cout << "k: " << k << endl;
                            convFlag = true;
                            break;
                        }
                    }
                    if (!convFlag)
                    {
                        cout << "FindClosestPointFromPolygonToMetaball2D: no convergence" << endl;
                        cout << "cp0: " << cp0.transpose() << endl;
                        cout << "initP: " << initP.transpose() << endl;
                        cout << "p0[i]: " << p0[i].transpose() << endl;
                        cout << "p0[id]: " << p0[id].transpose() << endl;
                        cout << "res: " << res << endl;
                        cout << "convk: " << convk << endl;
                        abort();
                    }
                    cp = cp0;
                    break;
                }
            }
            if (ci>maxC)
            {
                maxC = ci;
                maxID = id;
            }
            cb = ci;
            fb = fi;
        }
        if (!isLine)      cp = p0[maxID];
    }

    void FindMetaInitPointWithSurfacePoints(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, vector<Vector3d>& Sp0, vector<Vector3d>& Sp1, Vector3d& initP)
    {
        Vector3d cp0, cp1;
        FindClosestPointFromSurfacePoints(P1, K1, Sp0, cp0);
        FindClosestPointFromSurfacePoints(P0, K0, Sp1, cp1);
        initP = 0.5*(cp0+cp1);
    }

    void FindMetaClosestPoints(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, size_t nk, bool& conv, size_t& cvt)
    {
        Vector3d pk = initP;
        Vector3d dk = Vector3d::Zero();
        Vector3d dkn;

        Vector3d fp = Vector3d::Zero();
        Vector3d fq = Vector3d::Zero();
        // Newton's method
        conv = false;
        for (size_t k=0; k<=nk; ++k)
        {
            Vector3d f0 = Vector3d::Zero();
            Vector3d f1 = Vector3d::Zero();
            Matrix3d j0 = Matrix3d::Zero();
            Matrix3d j1 = Matrix3d::Zero();
            for (size_t i=0; i<P0.size(); ++i)
            {
                Vector3d a = pk-P0[i];
                double ri2inv = 1./a.squaredNorm();
                double ri4inv = ri2inv*ri2inv;
                double ri6inv = ri4inv*ri2inv;
                f0 -= 2.*K0(i)*ri4inv*a;
                j0(0,0) += 2.*K0(i)*(-ri4inv + 4.*a(0)*a(0)*ri6inv);
                j0(1,1) += 2.*K0(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
                j0(2,2) += 2.*K0(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
                j0(0,1) += 8.*K0(i)*a(0)*a(1)*ri6inv;
                j0(0,2) += 8.*K0(i)*a(0)*a(2)*ri6inv;
                j0(1,2) += 8.*K0(i)*a(1)*a(2)*ri6inv;
            }
            j0(1,0) = j0(0,1);
            j0(2,0) = j0(0,2);
            j0(2,1) = j0(1,2);
            for (size_t i=0; i<P1.size(); ++i)
            {
                Vector3d a = pk-P1[i];
                double ri2inv = 1./a.squaredNorm();
                double ri4inv = ri2inv*ri2inv;
                double ri6inv = ri4inv*ri2inv;
                f1 -= 2.*K1(i)*ri4inv*a;
                j1(0,0) += 2.*K1(i)*(-ri4inv + 4.*a(0)*a(0)*ri6inv);
                j1(1,1) += 2.*K1(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
                j1(2,2) += 2.*K1(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
                j1(0,1) += 8.*K1(i)*a(0)*a(1)*ri6inv;
                j1(0,2) += 8.*K1(i)*a(0)*a(2)*ri6inv;
                j1(1,2) += 8.*K1(i)*a(1)*a(2)*ri6inv;
            }
            j1(1,0) = j1(0,1);
            j1(2,0) = j1(0,2);
            j1(2,1) = j1(1,2);

            PartialPivLU<MatrixXd> lu(j0+j1);
            dkn = -lu.inverse()*(f0+f1);

            double dif = abs(dkn.norm()-dk.norm());

            if (dif>1.e-12)
            {
                dk = dkn;
                pk += dkn;
            }
            else
            {
                fp = f0;
                fq = f1;
                conv = true;
                cvt = k+1;
                break;
            }
        }

        double c0 = CalMetaC(P0, K0, pk);
        double c1 = CalMetaC(P1, K1, pk);

        Xc0 = pk+(1.-c0)/fp.squaredNorm()*fp;
        Xc1 = pk+(1.-c1)/fq.squaredNorm()*fq;

        finalP = pk;
    }
}