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

namespace Triangulation
{
    template <typename T>
    vector<VectorXi> TrianglizePolygon(const vector<T>& vertices) {
        auto isInside = [](const T& p, const T& a, const T& b, const T& c) {
            double area = 0.5 * (-b.y() * c.x() + a.y() * (-b.x() + c.x()) + a.x() * (b.y() - c.y()) + b.x() * c.y());
            double s = 1 / (2 * area) * (a.y() * c.x() - a.x() * c.y() + (c.y() - a.y()) * p.x() + (a.x() - c.x()) * p.y());
            double t = 1 / (2 * area) * (a.x() * b.y() - a.y() * b.x() + (a.y() - b.y()) * p.x() + (b.x() - a.x()) * p.y());
            return 0 <= s && s <= 1 && 0 <= t && t <= 1 && s + t <= 1;
        };

        auto isCCW = [](const T& a, const T& b, const T& c) {
            return (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()) > 0;
        };

        auto isEar = [&](const vector<size_t>& indices, size_t a, size_t b, size_t c) {
            if (!isCCW(vertices[indices[a]], vertices[indices[b]], vertices[indices[c]])) {
                return false;
            }
            for (size_t i = 0; i < indices.size(); ++i) {
                if (i != a && i != b && i != c) {
                    if (isInside(vertices[indices[i]], vertices[indices[a]], vertices[indices[b]], vertices[indices[c]])) {
                        return false;
                    }
                }
            }
            return true;
        };

        vector<VectorXi> faces;
        size_t n = vertices.size();

        if (n < 3) return faces;

        vector<size_t> indices(n);
        for (size_t i = 0; i < n; ++i) indices[i] = i;

        size_t current = 0;
        size_t remaining = n;

        while (remaining > 3) {
            size_t prev = (current + remaining - 1) % remaining;
            size_t next = (current + 1) % remaining;

            if (isEar(indices, prev, current, next)) {
                VectorXi face(3);
                face << indices[prev], indices[current], indices[next];
                faces.emplace_back(face);
                indices.erase(indices.begin() + current);
                --remaining;
                if (current == remaining) current = 0;
            } else {
                ++current;
                if (current == remaining) current = 0;
            }
        }
        VectorXi face(3);
        face << indices[0], indices[1], indices[2];
        faces.emplace_back(face);
        return faces;
    }

}