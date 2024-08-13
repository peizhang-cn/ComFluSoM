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

#ifndef DEM_CONTACT_CONVEX_METABALL_H
#define DEM_CONTACT_CONVEX_METABALL_H

inline void DEM::Metaball2Metaball(Vector3d xi, Vector3d xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, DEM_MAPS* maps, Vector3d& finalP, double& delta, Vector3d& n, Vector3d& cp)
{}

inline void DEM::Metaball2Wall(DEM_PARTICLE* Pi, DEM_PARTICLE* Pj, DEM_MAPS* maps, Vector3d& finalP, double& delta, Vector3d& n, Vector3d& cp)
{}

#endif