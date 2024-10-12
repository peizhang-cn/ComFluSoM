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

template<int SType, int D>
inline void MPM<SType, D>::CalNGN(MPM_PARTICLE* p0)
{
	double dx = Dx;
	Vector3d dxv (dx, dx, dx);
	Vector3d origin = Origin;
	Vector3d xp = p0->X;
	// Reset shape function (N) and gradient of shape function (GN)
	p0->ShapeInfo.resize(0);
	// Find min position of nodes which is infuenced by this particle
	Vector3i minx 	= Vector3i::Zero();
	Vector3i maxx 	= Vector3i::Zero();

	if constexpr (SType==1)
	{
		for (size_t d=0; d<D; ++d)
		{
			maxx(d) = (int) trunc((p0->X(d) + p0->PSize(d))/dx + 1.);
			minx(d) = (int) ceil((p0->X(d) - p0->PSize(d))/dx - 1.);
		}
	}
	else if constexpr (SType==2)
	{
		for (size_t d=0; d<D; ++d)
		{
			maxx(d) = (int) trunc(xp(d)/dx + 2.);
			minx(d) = (int) ceil(xp(d)/dx - 2.);
		}
	}

	// Find nodes within the influence range
	for (int i=minx(0); i<=maxx(0); ++i)
	for (int j=minx(1); j<=maxx(1); ++j)
	for (int k=minx(2); k<=maxx(2); ++k)
	{
		double n;
		Vector3d gn;
		Vector3d xn (i,j,k);
		xn = dx*xn + origin;
		// if constexpr (SType==GIMP && D==1)
		if constexpr (SType==1 && D==1)			MPM_ShapeFunction::GIMP1D(xp, xn, dxv, p0->PSize, n, gn);
		else if constexpr (SType==1 && D==2) 	MPM_ShapeFunction::GIMP2D(xp, xn, dxv, p0->PSize, n, gn);
		else if constexpr (SType==1 && D==3) 	MPM_ShapeFunction::GIMP3D(xp, xn, dxv, p0->PSize, n, gn);

		else if constexpr (SType==2)
		{
			cout << "\033[1;31mError: B-spline not supported yet! \033[0m\n";		
			exit(0);
		}

		if (n>0.)
		{
			int in, jn, kn;
			PeriodicNode(i,j,k,in,jn,kn);
			if (CheckIfInDomain(in,jn,kn))
			{
				size_t id = FindIDFromIJK(in,jn,kn);
				p0->ShapeInfo.push_back({id, n, gn});
				Ln[id]->Actived = true;
			}
		}
	}
}

template<int SType, int D>
inline void MPM<SType, D>::UpdateShapeFunction()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
	{
		MPM_PARTICLE* p0 = Lp[p];
		CalNGN(p0);
	}
}