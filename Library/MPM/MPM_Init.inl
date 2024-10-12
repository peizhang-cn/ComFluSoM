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
inline void MPM<SType, D>::Init(bool useFbar)
{
	cout << "================ Start init.  ================" << endl;
	Lp.resize(0);
	Ln.resize(0);
	Lcell.resize(0);
	// Add basic nodes
	for (size_t n=0; n<(Nx+1)*(Ny+1)*(Nz+1); ++n)
	{
    	size_t i, j, k;
    	FindIndex(n, i, j, k);
		Vector3d x (i, j, k);
		x *= Dx;
		Ln.push_back(new MPM_NODE(x));
    	Ln[Ln.size()-1]->ID = Ln.size()-1;
	}

	UseFbar = useFbar;
	cout << "=============== Finish init.  ================" << endl;
}

template<int SType, int D>
inline void MPM<SType, D>::InitCell()
{
	UseCell = true;
	for (size_t n=0; n<(Nx+1)*(Ny+1)*(Nz+1); ++n)
	{
    	size_t i, j, k;
    	FindIndex(n, i, j, k);
		Vector3d x (i+0.5, j+0.5, k+0.5);
		x *= Dx;
		Lcell.push_back(new MPM_CELL(x));
    	Lcell[Lcell.size()-1]->ID = Lcell.size()-1;
	}
}