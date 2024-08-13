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

#pragma once

class TRAJECTORY
{
public:
	TRAJECTORY(vector<Vector3d>& x, vector<Vector3d>& v, vector<Vector3d>& w);
	vector<Vector3d>	X;
	vector<Vector3d>	V;
	vector<Vector3d>	W;
};

inline TRAJECTORY::TRAJECTORY(vector<Vector3d>& x, vector<Vector3d>& v, vector<Vector3d>& w)
{
	if (x.size()!=v.size() || x.size()!=w.size() || v.size()!=w.size())
	{
		cout << "the size of TRAJECTORY x, v, w are not the same!" << endl;
		abort();
	}
	X = x;
	V = v;
	W = w;
}