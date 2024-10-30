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

namespace ConvexityCheck
{
    template<typename T>
    bool Polygon(const vector<T>& input_vertices) {
        static_assert(T::RowsAtCompileTime == 2 || T::RowsAtCompileTime == 3,
                    "This function only supports 2D or 3D vectors");

        int n = input_vertices.size();
        if (n < 3) return false;  // A polygon must have at least 3 vertices

        const bool is2D = (T::RowsAtCompileTime == 2);

        // Convert input to Vector3d
        vector<Vector3d> vertices;
        vertices.reserve(n);
        for (const auto& v : input_vertices) {
            if (is2D) {
                vertices.emplace_back(v.x(), v.y(), 0);
            } else {
                vertices.emplace_back(v);
            }
        }

        if (!is2D && n < 4) return false;  // A 3D polygon must have at least 4 vertices

        // Calculate the normal vector of the plane
        Vector3d normal = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]);

        if (!is2D) {
            // For 3D, check if all points lie on the same plane
            double d = normal.dot(vertices[0]);
            for (const auto& p : vertices) {
                if (abs(normal.dot(p) - d) > 1e-8) {
                    return false;  // Not all points are on the same plane
                }
            }
        }

        int sign = 0;
        for (int i = 0; i < n; ++i) {
            const Vector3d& p1 = vertices[i];
            const Vector3d& p2 = vertices[(i + 1) % n];
            const Vector3d& p3 = vertices[(i + 2) % n];

            Vector3d v1 = p2 - p1;
            Vector3d v2 = p3 - p2;
            double cross_product = v1.cross(v2).dot(normal);

            if (abs(cross_product) > 1e-8) {
                if (sign == 0) {
                    sign = (cross_product > 0) ? 1 : -1;
                } else if ((cross_product > 0 && sign < 0) || (cross_product < 0 && sign > 0)) {
                    return false;  // Found both positive and negative cross products
                }
            }
        }

        return true;  // All cross products have the same sign
    }
}