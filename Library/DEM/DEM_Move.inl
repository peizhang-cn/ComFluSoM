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

#ifndef DEM_MOVE_H
#define DEM_MOVE_H

inline void DEM::UpdatePositionGlobal(size_t t, double dt, bool resetHydro)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lp.size(); ++i)
	{
		DEM_PARTICLE* p0 = Lp[i];
		if (p0->isFixTrajectory)
		{
			size_t tid = p0->TrajID;
			p0->Xf = Ltraj[tid]->X[t];
			p0->Vf = Ltraj[tid]->V[t];
			p0->Wf = Ltraj[tid]->W[t];
		}
		if (p0->GroupID==-1)	p0->UpdatePosition(dt);
		p0->ZeroForceTorque(true, resetHydro);

		p0->CFMap.clear();
		p0->CFMap = std::move(p0->CFMapt);
		p0->CFMapt.clear();
		p0->NFMap.clear();
		p0->NFMap = std::move(p0->NFMapt);
		p0->NFMapt.clear();
		
		p0->MMap.clear();
		p0->MMap = std::move(p0->MMapt);
		p0->MMapt.clear();

		// for periodic BC
		for (size_t d=0; d<D; ++d)
		{
			if (Periodic[d])
			{
				if (p0->X(d)<Origin(d))			p0->X(d) += DomSize(d);
				else if (p0->X(d)>DomMax[d])	p0->X(d) -= DomSize(d);
			}
		}
		// p0->UpdateBox(D);

// discard REV function for now
/*		Vector3d xpRev = Tdi*p0->X;

		if (xpRev(0)>Lx)		xpRev(0) -= Lx;
		else if (xpRev(0)<0.)	xpRev(0) += Lx;
		if (xpRev(1)>Ly)		xpRev(1) -= Ly;
		else if (xpRev(1)<0.)	xpRev(1) += Ly;
		if (xpRev(2)>Lz)		xpRev(2) -= Lz;
		else if (xpRev(2)<0.)	xpRev(2) += Lz;

		p0->X = Td*xpRev;*/
	}

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lg.size(); ++i)
	{
		Lg[i]->UpdatePosition(dt);
		Lg[i]->P0->ZeroForceTorque(true, resetHydro);
	}
}

inline void DEM::UpdateVelocityGlobal(double dt)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lp.size(); ++i)
	{
		DEM_PARTICLE* p0 = Lp[i];
		if (p0->GroupID==-1)	p0->UpdateVelocity(dt);
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lg.size(); ++i)
	{
		Lg[i]->UpdateVelocity(dt);
	}
}

inline void DEM::ZeroForceTorqueGlobal(bool h, bool c)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->ZeroForceTorque(h, c);
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lg.size(); ++i)
	{
		Lg[i]->P0->ZeroForceTorque(h, c);
	}
}

#endif