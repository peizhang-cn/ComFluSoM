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

#ifndef DEM_PARTICLE_MOTION_H
#define DEM_PARTICLE_MOTION_H

inline void DEM_PARTICLE::UpdatePosition(double dt)
{
    if (isFixV)
    {
        V   = Vf;
        Avb.setZero();
    }
    
    if (isFixW)
    {
        W   = Wf;
        Awb.setZero();
    }

    X += dt*V + 0.5*dt*dt*Avb;
    //Update quaternion
    Quaterniond Aq, AAq, Qwb, Qawb;
    Qwb.w()     = 0.;
    Qwb.vec()   = 0.5*W;
    Qawb.w()    = 0.;
    Qawb.vec()  = 0.5*Awb;
    Quaterniond Qb = Q;
    // 5.43
    Aq   = Qb*Qwb;
    // Aq  = Qwb*Qb;
    // 5.44 and 5.46
    // AAq = Aq*Qwb + Qb*Qawb;
    // Q.coeffs() = Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*AAq.coeffs();

    Q.coeffs()   = Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Aq*Qwb).coeffs() + (Qb*Qawb).coeffs());
    // Q.coeffs()  = Qb.coeffs() + dt*Aq.coeffs() + 0.5*dt*dt*((Qwb*Aq).coeffs() + (Qawb*Qb).coeffs());
    Q.normalize();
    Qf = Q0*Q;  // final rotation (from object frame to lab frame)
    Qfi = Qf.inverse();
}

inline void DEM_PARTICLE::UpdateVelocity(double dt)
{
	Vector3d fht = Fh;
	Vector3d tht = Th;
    if (!isFixV)
    {
        // Acceleration at current time step
        Vector3d Av = (Fh + Fc + Fex)/M + G;
        // Update velocity
        V += 0.5*dt*(Avb+Av);
        // Store information for next update
        Avb = Av;
    }
    if (!isFixW)
    {
        // Total torque 
        Vector3d Tt = Th + Tc + Tex;
        //Update the angular velocity 5.51~53
        Vector3d Aw0 = I.asDiagonal().inverse()*Tt;
        // 5.45 and 54
        Vector3d w0 = W + 0.5*dt*(Awb + Aw0);
        //First order correction for angular velocity
        // 5.55-57
        Vector3d Aw1 = I.asDiagonal().inverse()*(-w0.cross(I.asDiagonal()*w0));
        // 5.58
        Vector3d w1 = w0 + 0.5*dt*Aw1;
        //Second order correction for angular velocity
        // 5.59-61
        Vector3d Aw2 = I.asDiagonal().inverse()*(-w1.cross(I.asDiagonal()*w1));
        // 5.62
        W   = w1 + 0.5*dt*Aw2;
        // Store information for next update
        // Awb = I.asDiagonal().inverse()*Tt;
		Awb(0) = (Tt(0) + (I(1)-I(2))*W(1)*W(2))/I(0);
		Awb(1) = (Tt(1) + (I(2)-I(0))*W(2)*W(0))/I(1);
		Awb(2) = (Tt(2) + (I(0)-I(1))*W(0)*W(1))/I(2);
		// if (Tt.norm()>1e-10)
		// {
		// 	cout << "W: " << W.transpose() << endl;
		// 	cout << "Tt: " << Tt.transpose() << endl;
		// 	cout << "===================" << endl;		    
		// }
    }
	if (isFixV)	Fh = fht;
	if (isFixW)	Th = tht;
}

inline void DEM_PARTICLE::ZeroForceTorque(bool h, bool c)
{
	//set zero hydro and contact force and torque after moved.
	if (h)
	{
		Fh.setZero();
		Th.setZero();
	}
	if (c)
	{
		Fc.setZero();
		Tc.setZero();
	}
}

inline void DEM_PARTICLE::UpdateBox(size_t D)
{
	//Update the range of wrapping box
	for (size_t d=0; d<D; ++d)
	{
		Max(d) = (int) (X(d)+BoxL(d));
		Min(d) = (int) (X(d)-BoxL(d));
	}
}

inline void DEM_PARTICLE::FixV(Vector3d& v)
{
	isFixV	= true;
	Vf	    = v;
}

inline void DEM_PARTICLE::FixW(Vector3d& w)
{
	isFixW	= true;
	Wf	= Q.inverse()._transformVector(w);
}

inline void DEM_PARTICLE::Fix()
{
	isFixed = true;
	Vector3d a = Vector3d::Zero();
	FixV(a);
	FixW(a);
}

inline void DEM_PARTICLE::UnFix()
{
	isFixed = false;
	isFixV = false;
	isFixW = false;
}

inline void DEM_PARTICLE::SetG(Vector3d& g)
{
    G       = g;
}

inline double DEM_PARTICLE::CalEnergy()
{
	double eng = 0.5*(M*V.squaredNorm()+I(0)*W(0)*W(0)+I(1)*W(1)*W(1)+I(2)*W(2)*W(2));
	return eng;
}

inline Vector3d DEM_PARTICLE::CalAngularMoment()
{
	Vector3d am = M*X.cross(V);
	return am;
}

inline void DEM_PARTICLE::AddReferencePoint(Vector3d x)
{
	Vector3d x0 = x-X;
	x0 = Q0.inverse()._transformVector(x0);
	P0.push_back(x0);
}

void inline DEM_PARTICLE::MooringBond(size_t ind, double l0, double kn, Vector3d xf)
{
	Vector3d point = GetP(ind);
	Vector3d xij = xf - point;
	double dis = xij.norm();
	Vector3d n = xij/dis;
	if (dis>l0)
	{
		Vector3d fn = kn*(dis-l0)*n;
		// cout << "fn: " << fn.transpose() << endl;
		Vector3d arm = point-X;
		Vector3d tc = Qfi._transformVector(arm).cross(Qfi._transformVector(fn));
		for (int d=0; d<3; ++d)
		{
			#pragma omp atomic
			Fc(d) += fn(d);
			#pragma omp atomic
			Tc(d) += tc(d);
		}
	}
}

inline Vector3d DEM_PARTICLE::GetP(size_t i)
{
	Vector3d p = Qf._transformVector(P0[i]) + X;
	// p = P0[i];
	return p;
}

inline Vector3d DEM_PARTICLE::Move2GlobalFrame(Vector3d xb)
{
	return (Qf._transformVector(xb) + X);
}

inline Vector3d DEM_PARTICLE::Move2BodyFrame(Vector3d x)
{
	return Qfi._transformVector(x - X);
}

#endif