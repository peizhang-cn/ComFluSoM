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

#include "../HEADER.h"

struct NEIGHBOR_LC
{
    size_t          ID;                     // ID of neighbor bin
    int             Periodic[3];            // Flag indicate if this is a periodic bin
};

struct BIN_LC
{
    size_t                              ID;         // ID in BINS
    Vector3d                            X;          // bin position
    vector<size_t>                      Lp;         // List of Particle ID that belongs to this bin
    vector<NEIGHBOR_LC>                 Ln;         // List of neighbor bin IDs at different levels bins
};

typedef pair<size_t, int[3]>    InteractPair;

class BINS_LC
{
public:
	BINS_LC(Vector3d origin, Vector3d l, Vector3i n, bool bx, bool by, bool bz);
    void SetParallel(size_t nproc);
    size_t FindBinID(Vector3d& x);                                                          // find bin ID by particle position x
    void FindBinIndex(size_t id, int& i, int& j, int& k);                                   // find i,j,k for bin from ID
    template <typename T>
    void AddParticleToBin(T& Par);                                                          // add particle to bin
    template <typename T>
    void UpdateBins(vector<T>& Lpar, bool& modified);                                       // update bins
    void FindInteractPair_P2P_Local(size_t p, size_t bid, vector<InteractPair>& Lc);        // find particle-particle contact pairs
    void FindInteractPair_P2W_Local(size_t p, vector<size_t> lb, vector<InteractPair>& Lc); // find particle-wall contact pairs

    size_t                              D;
    size_t                              Nproc;      // number of processers for parallel
    size_t                              Ncy;
    size_t                              Ncz;
    size_t                              Nb;         // total number of bins
    double                              MinL;       // min size of the bin
    Vector3d                            Origin;     // starting point of grid
    Vector3d                            Lx;         // domain size
    Vector3d                            Dx;         // size of each bin
    vector<int>                         NeiV;
    vector<BIN_LC*>                     Lb;         // list of bin
};

#include "BINS_UpdateBins.inl"
#include "BINS_FindInteractPair_P2W.inl"
#include "BINS_FindInteractPair_P2P.inl"

inline void BINS_LC::SetParallel(size_t nproc)
{
    Nproc = nproc;
}

inline size_t BINS_LC::FindBinID(Vector3d& x)
{
    Vector3d xr = x-Origin;
    size_t i = trunc(xr(0)/Dx(0));
    size_t j = trunc(xr(1)/Dx(1));
    size_t k = trunc(xr(2)/Dx(2));
    if (D==2)   k = 0;
    size_t id = i+j*Ncy+k*Ncz;
    return id;
}

inline void BINS_LC::FindBinIndex(size_t id, int& i, int& j, int& k)
{
    k = id/Ncz;
    j = (id%Ncz)/Ncy;
    i = (id%Ncz)%Ncy;
}

template <typename T>
inline void BINS_LC::AddParticleToBin(T& Par)
{
    size_t bid = FindBinID(Par->X);   // find bin id
    Lb[bid]->Lp.push_back(Par->ID);   // add to bin's Lp
}

BINS_LC::BINS_LC(Vector3d origin, Vector3d l, Vector3i n, bool bx, bool by, bool bz)
{
    Nproc = 1;
    Origin = origin;
    // lx,ly,lz: domain size, bx,by,bz: periodic flag for x,y,z-axis
    size_t nx = n(0);
    size_t ny = n(1);
    size_t nz = n(2);

    double lx = l(0);
    double ly = l(1);
    double lz = l(2);

    Lx = l;

    Dx(0) = lx/(double) nx;
    Dx(1) = ly/(double) ny;
    Dx(2) = lz/(double) nz;

    D = 3;
    if (lz==0.)
    {
        D = 2;
        // Dx(2) = 1.;
    }

    Ncz = nx*ny;                            // for FindBinID
    Ncy = nx;                               // for FindBinID
    if (D==2)           Nb = nx*ny;         // number of bin
    else if (D==3)      Nb = nx*ny*nz;      // number of bin

    if (D==2)           MinL = min(Dx(0), Dx(1));
    else if (D==3)      MinL = min(Dx(2), min(Dx(0), Dx(1)));

    vector<Vector3i>    Nei;                // Relative location (int) of neighbor bin
    if (D==3)
    {
        Nei  = {   { 0, 0, 0},
                   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
                   { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
                   { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
                   { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
                   { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
        NeiV.resize(27);
    }
    else if (D==2)
    {
        Nei  = {    { 0, 0, 0},
                    { 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
                    { 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
        NeiV.resize(9);
    }
    for (size_t i=0; i<NeiV.size(); ++i)    NeiV[i] = Nei[i](0)+Nei[i](1)*Ncy+Nei[i](2)*Ncz;

    Lb.resize(0);
    for (size_t b=0; b<Nb; ++b)             // generate bin system
    {
        BIN_LC* bin = new BIN_LC();
        bin->ID = b;
        bin->Lp.resize(0);
        bin->Ln.resize(0);
        Lb.push_back(bin);
    }

    int maxx[3] = {1,1,1};                  // for find bin neighbors
    if (D==2)   maxx[2] = 0;

    for (size_t b=0; b<Nb; ++b)             // add list of neighbor bins for each bin
    {
        int i, j, k;
        FindBinIndex(Lb[b]->ID, i, j, k);
        // find bin position
        Lb[b]->X(0) = ((double) i+0.5)*Dx(0);
        Lb[b]->X(1) = ((double) j+0.5)*Dx(1);
        Lb[b]->X(2) = ((double) k+0.5)*Dx(2);
        Lb[b]->X += origin;
        if (D==2)   Lb[b]->X(2) = 0.;

        for (int x=-maxx[0]; x<=maxx[0]; ++x)
        for (int y=-maxx[1]; y<=maxx[1]; ++y)
        for (int z=-maxx[2]; z<=maxx[2]; ++z)
        {
            int in = i+x; int jn = j+y; int kn = k+z;
            // avoid adding bin itself to neighbors
            bool itself = false;
            if (x==0 && y==0 && z==0)   itself = true;
            if (!itself)    // only for non-negative index
            {
                bool withInDom = true;
                int periodic[3] = {0, 0, 0};

                if (in == (int) nx)     // beyond maximum x
                {
                    if (bx) 
                    {
                        in = 0;
                        periodic[0] = 1;
                    }
                    else    withInDom = false;
                }
                else if (in == -1)      // beyond minimum x
                {
                    if (bx) 
                    {
                        in = nx-1;
                        periodic[0] = -1;
                    }
                    else    withInDom = false;
                }
                if (jn == (int) ny)     // beyond maximum y
                {
                    if (by) 
                    {
                        jn = 0;
                        periodic[1] = 1;
                    }
                    else    withInDom = false;
                }
                else if (jn == -1)      // beyond minimum y
                {
                    if (by) 
                    {
                        jn = ny-1;
                        periodic[1] = -1;
                    }
                    else    withInDom = false;
                }
                if (D==3)
                {
                    if (kn == (int) nz)     // beyond maximum z
                    {
                        if (bz) 
                        {
                            kn = 0;
                            periodic[2] = 1;
                        }
                        else    withInDom = false;
                    }
                    else if (kn == -1)      // beyond minimum z
                    {
                        if (bz) 
                        {
                            kn = nz-1;
                            periodic[2] = -1;
                        }
                        else    withInDom = false;
                    }
                }

                size_t idn = in+jn*Ncy+kn*Ncz;      // id for potential neighbor bin

                if (withInDom)
                {
                    // check if this bin already linked with current bin
                    bool exsit = false;
                    if (!exsit)                     // add if its a new link
                    {
                        NEIGHBOR_LC nei;
                        nei.ID = idn;
                        nei.Periodic[0] = periodic[0];
                        nei.Periodic[1] = periodic[1];
                        nei.Periodic[2] = periodic[2];
                        Lb[b]->Ln.push_back(nei);
                    }
                }
            }
        }
    }
}