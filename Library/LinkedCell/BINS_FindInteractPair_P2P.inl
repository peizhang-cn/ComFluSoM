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

inline void BINS_LC::FindInteractPair_P2P_Local(size_t p, size_t bid, vector<InteractPair>& Lc)
{
    BIN_LC* bin = Lb[bid];                      // current bin
    for (size_t i=0; i<bin->Lp.size(); ++i)     // for every particle in this bin
    {
        size_t q = bin->Lp[i];
        if (q>p)
        {
            InteractPair cp;
            cp.first = q;
            cp.second[0] = 0;
            cp.second[1] = 0;
            cp.second[2] = 0;
            Lc.push_back(cp);
        }
    }
    for (size_t n=0; n<bin->Ln.size(); ++n)  // loop over neighbor bins
    {
        size_t bn = bin->Ln[n].ID;           // neighbor bin ID
        BIN_LC* neighborBin = Lb[bn];           // neighbor bin
        size_t npn = neighborBin->Lp.size();    // numbers of particle in neighbor bin
        for (size_t j=0; j<npn; ++j)
        {
            size_t q = neighborBin->Lp[j];
            if (q>p)
            {
                InteractPair cp;
                cp.first = q;
                cp.second[0] = bin->Ln[n].Periodic[0];
                cp.second[1] = bin->Ln[n].Periodic[1];
                cp.second[2] = bin->Ln[n].Periodic[2];
                Lc.push_back(cp);
            }
        }
    }
}