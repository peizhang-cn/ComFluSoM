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

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

struct META_FUNC_FOR_MESH
{
	vector<Point_3>	MetaP;
	VectorXd		MetaK;
};

// warning global variable
META_FUNC_FOR_MESH* meta_para;

FT META_FUNC (Point_3 p) {
	FT c = -1.;
	for (size_t i=0; i<meta_para->MetaP.size(); ++i)
	{
		c += meta_para->MetaK(i)/squared_distance(p, meta_para->MetaP[i]);
	}
	return c;
}

// void CalSurfaceMesh(vector<Vector3d> metaP, VectorXd metaK, double dis, vector<Vector3d>& ps, vector<VectorXi>& faces)
// {
// 	meta_para = new META_FUNC_FOR_MESH();
// 	meta_para->MetaK = metaK;
// 	meta_para->MetaP.clear();

// 	Vector3d xc (0.,0.,0.);
// 	Vector3d maxL (-1.e30, -1.e30, -1.e30);
// 	Vector3d minL = -maxL;
// 	for (size_t i=0; i<metaP.size(); ++i)
// 	{
// 		xc += metaP[i];
// 		for (size_t j=0; j<3; ++j)
// 		{
// 			if (maxL(j)<metaP[i](j))	maxL(j) = metaP[i](j);
// 			if (minL(j)>metaP[i](j))	minL(j) = metaP[i](j);
// 		}
// 	}
// 	xc /= metaP.size();
// 	Vector3d length = maxL-minL;
// 	double maxLen = max(length(0), max(length(1), length(2)));
// 	double rescale = 1./maxLen;

// 	for (size_t i=0; i<metaP.size(); ++i)
// 	{
// 		metaP[i] -= xc;
// 		metaP[i] *= rescale;
// 		meta_para->MetaK(i) *= rescale*rescale;
// 	}

// 	Point_3 xc0 (0., 0., 0.);

// 	double rmax = 0.;
// 	for (size_t i=0; i<metaP.size(); ++i)
// 	{
// 		Point_3 mpi (metaP[i](0), metaP[i](1), metaP[i](2));
// 		meta_para->MetaP.push_back(mpi);
// 		double ri = metaP[i].norm();
// 		if (rmax<ri)	rmax = ri;
// 	}
// 	if (rmax==0.)	rmax = sqrt(metaK(0));
// 	// rmax *= 20.;
// 	cout << "rmax: " << rmax << endl;
// 	rmax *= 2.;
// 	cout << "rmax: " << rmax << endl;

// 	double dis_rescale = dis*rescale;

// 	Tr tr;            // 3D-Delaunay triangulation
// 	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
// 	// defining the surface
// 	Surface_3 surface(META_FUNC,             	// pointer to function
// 	                Sphere_3(xc0, rmax)); 		// bounding sphere

// 	// Note that "2." above is the *squared* radius of the bounding sphere!
// 	// defining meshing criteria
// 	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
// 	                                                 dis_rescale,  // radius bound
// 	                                                 dis_rescale); // distance bound	
// 	// meshing surface
// 	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
// 	Surface_mesh sm;
// 	CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

// 	ps.resize(0);
// 	for (Surface_mesh::Vertex_index vi : sm.vertices()) {
// 		Vector3d pi (sm.point(vi).x(), sm.point(vi).y(), sm.point(vi).z());
// 		pi /= rescale;
// 		pi += xc;
// 		ps.push_back(pi);
// 	}

// 	std::vector<int> indices;
// 	for (Surface_mesh::Face_index face_index : sm.faces()) {
// 	  CGAL::Vertex_around_face_circulator<Surface_mesh> vcirc(sm.halfedge(face_index), sm), done(vcirc);
// 	  do indices.push_back(*vcirc++); while (vcirc != done);
// 	}

// 	faces.resize(0);
// 	for (size_t i=0; i<(indices.size()+1)/3; ++i)
// 	{
// 		VectorXi facei;
// 		facei.resize(3);
// 		facei(0) = indices[3*i];
// 		facei(1) = indices[3*i+1];
// 		facei(2) = indices[3*i+2];
// 		faces.push_back(facei);
// 	}

// 	cout << "ps: " << ps.size() << endl;
// 	cout << "faces: " << faces.size() << endl;

// 	delete	meta_para;
// }

void CalSurfaceMesh(vector<Vector3d> metaP, VectorXd metaK, double dis, vector<Vector3d>& ps, vector<VectorXi>& faces)
{
	meta_para = new META_FUNC_FOR_MESH();
	meta_para->MetaK = metaK;
	meta_para->MetaP.clear();

	Vector3d xc (0.,0.,0.);
	for (size_t i=0; i<metaP.size(); ++i)	xc += metaP[i];
	xc /= metaP.size();

	Point_3 xc0 (xc(0), xc(1), xc(2));

	double rmax = 0.;
	for (size_t i=0; i<metaP.size(); ++i)
	{
		Point_3 mpi (metaP[i](0), metaP[i](1), metaP[i](2));
		meta_para->MetaP.push_back(mpi);
		double ri = (metaP[i]-xc).norm();
		if (rmax<ri)	rmax = ri;
	}
	if (rmax==0.)	rmax = sqrt(metaK(0));
	rmax *= 10.;

	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
	// defining the surface
	Surface_3 surface(META_FUNC,             	// pointer to function
	                Sphere_3(xc0, rmax)); 		// bounding sphere

	// Note that "2." above is the *squared* radius of the bounding sphere!
	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
	                                                 dis,  // radius bound
	                                                 dis); // distance bound	
	// meshing surface
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	Surface_mesh sm;
	CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

	ps.resize(0);
	for (Surface_mesh::Vertex_index vi : sm.vertices()) {
		Vector3d pi (sm.point(vi).x(), sm.point(vi).y(), sm.point(vi).z());
		ps.push_back(pi);
	}

	std::vector<int> indices;
	for (Surface_mesh::Face_index face_index : sm.faces()) {
	  CGAL::Vertex_around_face_circulator<Surface_mesh> vcirc(sm.halfedge(face_index), sm), done(vcirc);
	  do indices.push_back(*vcirc++); while (vcirc != done);
	}

	faces.resize(0);
	for (size_t i=0; i<(indices.size()+1)/3; ++i)
	{
		VectorXi facei;
		facei.resize(3);
		facei(0) = indices[3*i];
		facei(1) = indices[3*i+1];
		facei(2) = indices[3*i+2];
		faces.push_back(facei);
	}

	cout << "ps: " << ps.size() << endl;
	cout << "faces: " << faces.size() << endl;

	delete	meta_para;
}