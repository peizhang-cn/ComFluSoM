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

#ifndef DEM_SOLVE_H
#define DEM_SOLVE_H

inline void DEM::SolveOneStep(size_t t, size_t ts, vector<double>& times)
{
	bool firststep = false;
	// bool save = false;
	if (t%ts == 0)
	{
		// save = true;
		WriteFileH5(t);
	}
	if (t==0)
	{
		for (size_t p=0; p<Lp.size(); ++p)
		{
			Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
			Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
		}
		firststep = true;
	}

	double time0 = 0.;
	auto t_start = std::chrono::system_clock::now();
	UpdatePositionGlobal(t, Dt, true);
	auto t_end = std::chrono::system_clock::now();
	time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
	// if (show)
	{
		// cout << "UpdatePositionGlobal time= " << time0 << endl;
		times.push_back(time0);
		time0 = 0.;
	}
	// cout << "Lp[0]->X: " << Lp[0]->X.transpose() << endl;
	// cout << "Lp[1]->X: " << Lp[1]->X.transpose() << endl;

	// cout << "Lp[0]->V: " << Lp[6]->V.transpose() << endl;
	// cout << "Lp[1]->V: " << Lp[1]->V.transpose() << endl;

	// cout << "Lp[0]->W: " << Lp[0]->W.transpose() << endl;
	// cout << "Lp[1]->W: " << Lp[1]->W.transpose() << endl;

	t_start = std::chrono::system_clock::now();
	LinkedCell(firststep);
	t_end = std::chrono::system_clock::now();
	time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
	// if (show)
	{
		// cout << "LinkedCell time= " << time0 << endl;
		times.push_back(time0);
		time0 = 0.;
	}

	// cout << "Lc: " << Lc.size() << endl;
	t_start = std::chrono::system_clock::now();
	Contact(false, 0, false);
	t_end = std::chrono::system_clock::now();
	time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
	// if (show)
	{
		// cout << "Contact time= " << time0 << endl;
		times.push_back(time0);
		time0 = 0.;
	}

	t_start = std::chrono::system_clock::now();
	UpdateVelocityGlobal(Dt);
	t_end = std::chrono::system_clock::now();
	time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
	// if (show)
	{
		// cout << "UpdateVelocityGlobal time= " << time0 << endl;
		times.push_back(time0);
		time0 = 0.;
		// cout << "===========================" << endl;
	}
}

// inline void DEM::Solve(int tt, int ts, bool writefile)
// {
// 	if (!BinExist)
// 	{
// 		cout << "\033[1;31mError: Bin system does not exit! InitBinSystem first before Solve! \033[0m\n";		
// 		exit(0);			
// 	}
// 	double time0;
// 	time0 = 0.;
// 	for (int t=0; t<tt; ++t)
// 	{
// 		bool show = false;
// 		bool save = false;
// 		bool firststep = false;
// 		if (t%ts == 0)
// 		{
// 			show = true;
// 			save = true;
// 			cout << "Time Step ============ " << t << endl;
// 		}
// 		if (t==0)
// 		{
// 			for (size_t p=0; p<Lp.size(); ++p)
// 			{
// 				Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
// 				Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
// 			}
// 			firststep = true;
// 		}		
// 		auto t_start = std::chrono::system_clock::now();
// 		UpdatePositionGlobal(t, Dt, true);
// 		auto t_end = std::chrono::system_clock::now();
// 		time0 = 0.;
// 		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		if (show)
// 		{
// 			cout << "UpdatePositionGlobal time= " << time0 << endl;
// 			time0 = 0.;
// 		}

// 		size_t chunk_size = 8;
// 		t_start = std::chrono::system_clock::now();
// 		LinkedCell(firststep, chunk_size);
// 		t_end = std::chrono::system_clock::now();
// 		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		if (show)
// 		{
// 			cout << "LinkedCell time= " << time0 << endl;
// 			time0 = 0.;
// 		}

// 		t_start = std::chrono::system_clock::now();
// 		Contact(false, 0, show);
// 		t_end = std::chrono::system_clock::now();
// 		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		if (show)
// 		{
// 			cout << "Contact time= " << time0 << endl;
// 			time0 = 0.;
// 			cout << "NmetaContact: " << NmetaContact << endl;
// 			cout << "Niteration: " << Niteration << endl;
// 			cout << "Niteration/NmetaContact: " << (double) Niteration/ (double) NmetaContact << endl;
// 		}
// 		NmetaContact = 0;
// 		Niteration = 0;

// 		t_start = std::chrono::system_clock::now();
// 		UpdateVelocityGlobal(Dt);
// 		t_end = std::chrono::system_clock::now();
// 		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		if (show)
// 		{
// 			cout << "UpdateVelocityGlobal time= " << time0 << endl;
// 			time0 = 0.;
// 		}
// 		if (writefile && save)	WriteFileH5(t);
// 	}
// }

#endif