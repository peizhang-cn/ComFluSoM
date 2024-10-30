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
inline void MPM<SType, D>::SolveMUSL(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (ActiveShift)
		{
			for (size_t d=0; d<D; ++d)	Shift(d) = GetUniformD1();		
			// for (size_t d=0; d<D; ++d)	Shift(d) = (t%2)*0.5;		
		}
		bool show = false;
		if (t%ts==0)	show = true;
		if (show) 	cout << "Time Step = " << t << endl;
		if (t%ts == 0)
		{
			// WriteFileH5(t);
			WriteFileH5Particle("",t);
		}

		auto t_start = std::chrono::system_clock::now();
		UpdateTopographicForce(1.e6, 0.);
		auto t_end = std::chrono::system_clock::now();
		if (show)	cout << "UpdateTopologyFc= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		// UpdateShapeFunction();
		ParticleToNode();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "ParticleToNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		CalVOnNode();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalVOnNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		NodeToParticle();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "NodeToParticle= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		ResetNodes();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "ResetNodes= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (show) 	cout << "===========================" << endl;
	}
}