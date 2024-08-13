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

#ifndef DEM_CONTACT_DETECTION_H
#define DEM_CONTACT_DETECTION_H

inline void DEM::LinkedCell(bool firststep)
{
	bool modified = false;
	Bins->UpdateBins(Lp, modified);

	if (modified || firststep)
	{
		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t p=0; p<Lp.size(); ++p)
		{
			DEM_PARTICLE* par = Lp[p];
			if (!par->isBig)
			{
				par->Lc.clear();
				size_t bid = Bins->FindBinID(par->X); 	// find bin id
				Bins->FindInteractPair_P2P_Local(p, bid, par->Lc);
			}
		}
	}

	double maxDis = Bins->Dx.norm();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lw.size(); ++i)
	{
        size_t p = Lw[i];
        DEM_PARTICLE* par = Lp[p];
        if (!par->isFixed || firststep)
        {
            par->Lb.clear();
            for (size_t b=0; b<Bins->Nb; ++b)
            {
                Vector3d xb = Bins->Lb[b]->X;
                double dis = abs(par->GetDistance(xb));
                if (dis<maxDis) par->Lb.push_back(b);
            }
        }
	}

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Lw.size(); ++i)
	{
        size_t p = Lw[i];
        DEM_PARTICLE* par = Lp[p];
        if (!par->isFixed || firststep || modified)
        {
			// cout << "update lc" << endl;
	        par->Lc.clear();
	        Bins->FindInteractPair_P2W_Local(p, par->Lb, par->Lc);
			// cout << "par->Lc: " << par->Lc.size() << endl;
        }
	}
}

#endif