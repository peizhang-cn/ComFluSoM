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

template<int SType, int D>
inline void MPM<SType, D>::SetTopography(vector<Vector3d> vertices, vector<VectorXi> faces)
{
	UseTopology = true;
	Topology->Vertices = vertices;
	Topology->Faces = faces;

	if (D==2)
	{
		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t n=0; n<Ln.size(); ++n)
		{
			Vector3d xc = Ln[n]->X;
			double sdf = SignedDistance::Polygon2D(xc, vertices, faces);
			Ln[n]->SDF = sdf;
		}
	}
	else if (D==3)
	{
		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t n=0; n<Ln.size(); ++n)
		{
			Vector3d xc = Ln[n]->X;
			double sdf = SignedDistance::Polyhedron3D(xc, vertices, faces);
			Ln[n]->SDF = sdf;
		}
	}
}

template<int SType, int D>
inline double MPM<SType, D>::CalParticleSDF(MPM_PARTICLE* p0)
{
	double sdf = 0.;
	CalNGN(p0);
	for (size_t l=0; l<p0->ShapeInfo.size(); ++l)
	{
		size_t id = p0->ShapeInfo[l].ID;
		double n = p0->ShapeInfo[l].N;
		sdf += n*Ln[id]->SDF;
	}
	return sdf;
}

template<int SType, int D>
inline void MPM<SType, D>::UpdateParticleTopographicForce(MPM_PARTICLE* p0, double kn, double gn)
{
	double sdf = 0.;
	Vector3d gsdf = Vector3d::Zero();
	for (size_t l=0; l<p0->ShapeInfo.size(); ++l)
	{
		size_t id = p0->ShapeInfo[l].ID;
		double n = p0->ShapeInfo[l].N;
		Vector3d gn = p0->ShapeInfo[l].Gn;
		sdf += n*Ln[id]->SDF;
		gsdf += Ln[id]->SDF*gn;
	}
	double delta = p0->PSize0(0)-sdf;
	if (delta>0.)
	{
		gsdf.normalized();
		p0->Fc += kn*delta*gsdf;
		// p0->M = 0.;
	}
}

template<int SType, int D>
inline void MPM<SType, D>::UpdateTopographicForce(double kn, double gn)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateParticleTopographicForce(Lp[p], kn, gn);
	}
}

template<int SType, int D>
void MPM<SType, D>::SetNonSlippingBC(size_t n)
{
	Ln[n]->BCTypes.push_back(1);
}

template<int SType, int D>
void MPM<SType, D>::SetNonSlippingBC(size_t i, size_t j, size_t k)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(1);
}

template<int SType, int D>
void MPM<SType, D>::SetSlippingBC(size_t n, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(2);
	Ln[n]->Norms.push_back(norm);
}

template<int SType, int D>
void MPM<SType, D>::SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	// cout << "n= " << n << endl;
	// cout << "Nnode= " << Nnode << endl;
	Ln[n]->BCTypes.push_back(2);
	// cout << "push_back 1 " << endl;
	Ln[n]->Norms.push_back(norm);
	// cout << "push_back 2 " << endl;
}

template<int SType, int D>
void MPM<SType, D>::SetFrictionBC(size_t n, double mu, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}

template<int SType, int D>
void MPM<SType, D>::SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}