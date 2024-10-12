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
void MPM<SType, D>::SetNonSlippingBC(size_t n)
{
	Ln[n]->BCTypes.push_back(1);
}

template<int SType, int D>
void MPM<SType, D>::SetNonSlippingBC(size_t i, size_t j, size_t k)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(1);
}

template<int SType, int D>
void MPM<SType, D>::SetSlippingBC(size_t n, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(2);
	Ln[n]->Norms.push_back(norm);
}

template<int SType, int D>
void MPM<SType, D>::SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	// cout << "n= " << n << endl;
	// cout << "Nnode= " << Nnode << endl;
	Ln[n]->BCTypes.push_back(2);
	// cout << "push_back 1 " << endl;
	Ln[n]->Norms.push_back(norm);
	// cout << "push_back 2 " << endl;
}

template<int SType, int D>
void MPM<SType, D>::SetFrictionBC(size_t n, double mu, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}

template<int SType, int D>
void MPM<SType, D>::SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}