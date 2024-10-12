/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2021 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License fo`r more details.                        *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/*void MPM::ParticleToCell()
{
   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t n=0; n<Lcell.size(); ++n)
   {
   	Lcell[n]->Reset();
   }

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		MPM_PARTICLE* p0 = Lp[p];
		int i = ceil(p0->X(0))-1;
		int j = ceil(p0->X(1))-1;
		int k = ceil(p0->X(2))-1;
		// Find id of current cell
		size_t ii = (i+Nx+1)%(Nx+1);
		size_t jj = (j+Ny+1)%(Ny+1);
		size_t kk = (k+Nz+1)%(Nz+1);
		size_t cid = ii+jj*Ncy+kk*Ncz;
		p0->CID = cid;

		Vector3d xcp = Lcell[cid]->X-p0->X;
		if (abs(xcp(0))>0.5 || abs(xcp(1))>0.5)
		{
			cout << "Lcell[cid]->X: " << Lcell[cid]->X.transpose() << endl;
			cout << "p0->X: " << p0->X.transpose() << endl;
			abort();
		}

		#pragma omp atomic
		Lcell[cid]->Vol0 += p0->Vol0;
		#pragma omp atomic
		Lcell[cid]->Vol += p0->Vol0*p0->Td.determinant();

		// #pragma omp atomic
		// Lcell[cid]->Np++;

		// for (size_t d=0; d<D; ++d)
		// for (size_t c=0; c<D; ++c)
		// {
		// 	#pragma omp atomic
		// 	Lcell[cid]->Stress(d,c) += p0->Stress(d,c);
		// }		
	}

   #pragma omp parallel for schedule(static) num_threads(Nproc)
   for (size_t n=0; n<Lcell.size(); ++n)
   {
   	Lcell[n]->J = Lcell[n]->Vol/Lcell[n]->Vol0;
   	// cout << "j: " << Lcell[n]->J << endl;
   	// Lcell[n]->Stress /= (double) Lcell[n]->Np;
   }
}*/