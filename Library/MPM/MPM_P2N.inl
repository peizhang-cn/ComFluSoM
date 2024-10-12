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
inline void MPM<SType, D>::ParticleToNode()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	MPM_PARTICLE* p0 = Lp[p];
		CalNGN(p0);
    	Matrix3d vsp = -p0->Vol*p0->Stress;
		Vector3d fex = p0->M*p0->B + p0->Fh + p0->Fc;
		for (size_t l=0; l<p0->ShapeInfo.size(); ++l)
		{
			// Grid id
			size_t id = p0->ShapeInfo[l].ID;
			// weight
			double 		n 	= p0->ShapeInfo[l].N;
			Vector3d 	gn 	= p0->ShapeInfo[l].Gn;
			Vector3d 	df 	= n*fex + vsp*gn;
			// weigthed mass contribution
			double nm = n*p0->M;
			Vector3d ve = p0->V;
			// Vector3d ve = p0->V-Lp[p]->L*(Lp[p]->X-Ln[id]->X);
			#pragma omp atomic
			Ln[id]->M += nm;

			#pragma omp atomic
			Ln[id]->Vol += n*p0->Vol;
			#pragma omp atomic
			Ln[id]->JbarDeltaJ += n*p0->JbarDeltaJ*p0->Vol;

			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->Mv(d) += nm*ve(d);
				#pragma omp atomic
				Ln[id]->F(d) += df(d);
			}
		}
    }

	vector<vector <size_t>> lan(Nproc);
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Ln.size(); ++n)
	{
		if (Ln[n]->Actived)
		{
			auto id = omp_get_thread_num();
			lan[id].push_back(n);
		}
	}

	for (size_t n=0; n<Nproc; ++n)
	{
		LAn.insert( LAn.end(), lan[n].begin(), lan[n].end() );
	}
}