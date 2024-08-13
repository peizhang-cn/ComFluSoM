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

#ifndef __TRUABGLIZAPOLYGON_H__
#define __TRUABGLIZAPOLYGON_H__

namespace TrianglizePolygon
{
    template <typename T>
    void TrianglizePolygon(vector<T>& P, vector<VectorXi>& F)
    {
        vector<pair<size_t, T>> pt;
        for (size_t i = 0; i < P.size(); ++i)
        {
            pt.push_back(make_pair(i,P[i]));
        }

        while (pt.size()>3)
        {
            int rmID = 0;
            for (int i = 0; i < pt.size(); ++i)
            {
                int ip = (i-1+pt.size())%pt.size();
                int in = (i+1+pt.size())%pt.size();
                T vec1 = pt[ip].second-pt[i].second;
                T vec2 = pt[in].second-pt[i].second;
                double cross = vec1(0) * vec2(1) - vec1(1) * vec2(0);
                if (cross<=0.)
                {
                    VectorXi face(3);
                    face << pt[ip].first, pt[i].first, pt[in].first;
                    F.push_back(face);
                    rmID = i;
                    break;
                }
            }
            pt.erase(pt.begin() + rmID);
        }
        VectorXi face(3);
        face << pt[0].first, pt[1].first, pt[2].first;
        F.push_back(face);
    }
}

#endif