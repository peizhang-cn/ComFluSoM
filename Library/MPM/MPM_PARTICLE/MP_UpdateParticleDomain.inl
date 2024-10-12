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

inline void MPM_PARTICLE::UpdatePSize(size_t flag)
{
	switch (flag)
	{
	case 1:
		{
			PSize(0) =  PSize0(0)*Td(0,0);
			PSize(1) =  PSize0(1)*Td(1,1);
			PSize(2) =  PSize0(2)*Td(2,2);
		}
	case 2:
		{
			PSize(0) =  PSize0(0)*sqrt(Td(0,0)*Td(0,0) + Td(1,0)*Td(1,0) + Td(2,0)*Td(2,0));
			PSize(1) =  PSize0(1)*sqrt(Td(0,1)*Td(0,1) + Td(1,1)*Td(1,1) + Td(2,1)*Td(2,1));
			PSize(2) =  PSize0(2)*sqrt(Td(0,2)*Td(0,2) + Td(1,2)*Td(1,2) + Td(2,2)*Td(2,2));		
		}
	}
}