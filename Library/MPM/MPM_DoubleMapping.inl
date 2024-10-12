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

/*void MPM::NodeToParticleDoubleMapping()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateParticle(Lp[p]);
	}
	CalVOnNodeDoubleMapping();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateStress(Lp[p]);
	}		
}*/

/*void MPM::CalVOnNodeDoubleMapping()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<LAn.size(); ++c)
	{
		Ln[c]->V = Vector3d::Zero();
	}
	// Double map velocity from particles to nodes to aviod small mass problem
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 	n = Lp[p]->LnN[l];
			// weigthed mass contribution
			double nm = n*Lp[p]->M;
			// #pragma omp atomic
			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->V(d) += nm*Lp[p]->V(d)/Ln[id]->M;
			}
		}
    }
	// #pragma omp parallel for schedule(static) num_threads(1)
	// for (size_t c=0; c<LAn.size(); ++c)
	// {
	// 	Ln[c]->V = Ln[c]->Mv/Ln[c]->M;
	// }
}*/
