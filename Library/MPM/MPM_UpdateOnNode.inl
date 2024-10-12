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
void MPM<SType, D>::CalVOnNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Ln[id]->Stress /= Ln[id]->M;

		Ln[id]->F  *= (1.-Dc*MPM_ShapeFunction::Sign(Ln[id]->Mv.dot(Ln[id]->F)));
		Ln[id]->Mv += Ln[id]->F*Dt;

		if (Ln[id]->BCTypes.size()>0)
		{
			for (size_t i=0; i<Ln[id]->BCTypes.size(); ++i)
			{
				if (Ln[id]->BCTypes[i]==1)			Ln[id]->NonSlippingBC();
				else if (Ln[id]->BCTypes[i]==2)		Ln[id]->SlippingBC(Ln[id]->Norms[i]);
				else if (Ln[id]->BCTypes[i]==3)		Ln[id]->FrictionBC(Dt, Ln[id]->Norms[i]);
			}
		}
		Ln[id]->V = Ln[id]->Mv/Ln[id]->M;
		Ln[id]->F /= Ln[id]->M;

		Ln[id]->JbarDeltaJ /= Ln[id]->Vol;
	}
}