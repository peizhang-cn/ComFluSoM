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

template <typename T>
inline void BINS_LC::UpdateBins(vector<T>& Lpar, bool& modified)
{
    modified = false;
    // Loop for all bins and find particles that need to update their bin ID
    vector< vector<pair<size_t,size_t>> > updateLpProc(Nproc);      // list of particles that need to update it's bin id
    for (size_t c=0; c<Nproc; ++c)
    {
        updateLpProc[c].clear();
    }

    #pragma omp parallel for schedule(static, 1) num_threads(Nproc)
    for (size_t b=0; b<Nb; ++b)   // loop all bins
    {
        BIN_LC* bin = Lb[b];      // current bin
        for (size_t i=0; i<bin->Lp.size(); ++i)
        {
            size_t j = bin->Lp.size()-1-i;  // loop from the last particle in the bin
            size_t p = bin->Lp[j];          // id for particle belongs to this bin

            size_t bid = FindBinID(Lpar[p]->X); // find bin id
            if (bid != b)
            {
                modified = true;
                swap(bin->Lp[j], bin->Lp.back());       // move to end of Lp
                bin->Lp.pop_back();                     // remove particle from bin
                pair<size_t, size_t> bid_p;             // a pair that first is the bin ID and second is ID of particle that need to be added.
                bid_p.first = bid;
                bid_p.second = p;
                int tid = omp_get_thread_num();      // thread id
                updateLpProc[tid].push_back(bid_p);
            }
        }
    }
    // Add particles to their bin
    // notice that by assuming there is only very few particles need to be move, therefore, this part is not parallelized
    for (size_t c=0; c<Nproc; ++c)
    for (size_t i=0; i<updateLpProc[c].size(); ++i)
    {
        size_t bid = updateLpProc[c][i].first;
        size_t p = updateLpProc[c][i].second;
        Lb[bid]->Lp.push_back(p);         // add to bin's Lp
    }
}