/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#include "../HEADER.h"
#include <SHAPE.h>
#include <MPM_PARTICLE.h>
#include <MPM_NODE.h>

class MPM
{
public:
	MPM();
	~MPM();
	MPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx);
	void Init();
	void UpdateLn(MPM_PARTICLE* p0);
	void CalNGN(MPM_PARTICLE* p0);
	void CalNGN_MLS(MPM_PARTICLE* p0);
	void UpdateLAn();
	void ParticleToNode();
	void ParticleToNodeMLS();
	void CalFOnNode(bool firstStep);
	void NodeToParticle();
	void CalStressOnParticleElastic();
	void CalStressOnParticleMohrCoulomb();
	void CalStressOnParticleNewtonian();
	void CalVGradLocal(int p);
	void CalPSizeCP(int p);
	void CalPSizeR(int p);
	void CalVOnNode();
	void CalVOnNodeDoubleMapping();
	void SetNonSlippingBC(size_t n);
	void SetNonSlippingBC(size_t i, size_t j, size_t k);
	void SetSlippingBC(size_t n, Vector3d& norm);
	void SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm);
	void SetFrictionBC(size_t n, double mu, Vector3d& norm);
	void SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm);
	void SolveMUSL(int tt, int ts);
	void SolveUSF(int tt, int ts);
	void SolveUSA(int tt, int ts);
	void AddNode(size_t level, Vector3d& x);
	// void AddParticle(int tag, Vector3d& x, double m, double young, double poisson);
	void AddParticle(int tag, Vector3d& x, double m);
	void DeleteParticles();
	// void AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson);
	void AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m);
	void WriteFileH5(int n);
	void FindIndex(size_t n, size_t& i, size_t& j, size_t& k);

	double 		(*N)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	Vector3d 	(*GN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp);
	void 		(*NGN)(Vector3d& x, Vector3d& xc, Vector3d& l, Vector3d& lp, double& n, Vector3d& gn);

	vector <size_t>					LAn;													// List of actived nodes
	vector <MPM_PARTICLE*>			Lp;														// List of all MPM particles
	vector <MPM_PARTICLE*>			Lbp;													// List of boundary MPM particles
	vector <MPM_NODE*>				Ln;														// List of all MPM nodes

	bool							Periodic[3];
	bool 							MLSv;

    size_t 							Nx;														// Domain size
    size_t 							Ny;
    size_t 							Nz;
    size_t 							Ncz;
    size_t 							Ncy;
    size_t 							Nnode;													// Total number of nodes

    size_t 							Nproc;
    size_t 							D;														// Dimension	
    size_t 							Ntype;													// Type of shape function 0 for Linear, 1 for Quadratic and 2 for Cubic 3 for GIMP

    double 							Nrange;													// Influence range of shape function
    double 							Dt;														// Time step
    double 							Dc;														// Damping coefficient
    double 							Cs;														// Speed of sound
    Vector3d						Dx;														// Space step
};

MPM::MPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx)
{
	Nproc	= 1;

	Ntype 	= ntype;

	Nx 		= nx;
	Ny 		= ny;
	Nz 		= nz;
	Dx 		= dx;
	D 		= 3;
	Dt 		= 1.;
	Dc 		= 0.;
	Cs 		= 0.;
	Periodic[0] = false;
	Periodic[1] = false;
	Periodic[2] = false;

	MLSv = false;

	Ncz = (Nx+1)*(Ny+1);
	Ncy = (Nx+1);

	Nnode = (Nx+1)*(Ny+1)*(Nz+1);

	if (Nz==0)
	{
		D = 2;
		if (Ny==0)	D = 1;
	}
	// Linear
	if 		(Ntype == 0)
	{
		if (D==1) 			NGN =& LS1D;
		else if (D==2)		NGN =& LS2D;
		else if (D==3)		NGN =& LS3D;
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.;
		cout << "Using Linear shape function." << endl;
	}
	// // Quadratic
	// else if (Ntype == 1)
	// {
	// 	if (D==1)
	// 	{
	// 		N  		=& ShapeQ1D;
	// 		GN 		=& GradShapeQ1D;			
	// 	}
	// 	else if (D==2)
	// 	{
	// 		N  		=& ShapeQ2D;
	// 		GN 		=& GradShapeL2D;
	// 	}
	// 	else if (D==3)
	// 	{
	// 		N  		=& ShapeQ3D;
	// 		GN 		=& GradShapeQ3D;
	// 	}
	// 	else
	// 	{
	// 		cout << "Dimension is higher than 3." << endl;
	// 		abort();
	// 	}
	// 	Nrange 	= 1.5;
	// 	cout << "Using Quadratic shape function." << endl;
	// }
	// // Cubic
	// else if (Ntype == 2)
	// {
	// 	if (D==1)
	// 	{
	// 		N  		=& ShapeC1D;
	// 		GN 		=& GradShapeC1D;			
	// 	}
	// 	else if (D==2)
	// 	{
	// 		N  		=& ShapeC2D;
	// 		GN 		=& GradShapeC2D;
	// 	}
	// 	else if (D==3)
	// 	{
	// 		N  		=& ShapeC3D;
	// 		GN 		=& GradShapeC3D;
	// 	}
	// 	else
	// 	{
	// 		cout << "Dimension is higher than 3." << endl;
	// 		abort();
	// 	}
	// 	Nrange 	= 2.;
	// 	cout << "Using Cubic shape function." << endl;
	// }
	// GIMP
	else if (Ntype == 3)
	{
		if (D==1)			NGN =& GIMP1D;
		else if (D==2) 		NGN =& GIMP2D;
		else if (D==3)		NGN =& GIMP3D;
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.;
		cout << "Using GIMP shape function." << endl;
	}
	else
	{
		cout << "Undefined shape function type. Retry 0 for Linear, 1 for Quadratic and 2 for Cubic." << endl;
		abort();
	}
}

void MPM::Init()
{
	cout << "================ Start init.  ================" << endl;
	Lp.resize(0);
	Ln.resize(0);
	// Add basic nodes
	for (size_t n=0; n<(Nx+1)*(Ny+1)*(Nz+1); ++n)
	{
    	size_t i, j, k;
    	FindIndex(n, i, j, k);
		Vector3d x (i, j, k);
		Ln.push_back(new MPM_NODE(0,x));
    	Ln[Ln.size()-1]->ID = Ln.size()-1;
	}
	cout << "=============== Finish init.  ================" << endl;
}
// Find index for grid
inline void MPM::FindIndex(size_t n, size_t& i, size_t& j, size_t& k)
{
	k = n/Ncz;
	j = (n%Ncz)/Ncy;
	i = (n%Ncz)%Ncy;
}

void MPM::CalNGN(MPM_PARTICLE* p0)
{
	// Reset shape function (N) and gradient of shape function (GN)
	p0->Lni.resize(0);
	p0->LnN.resize(0);
	p0->LnGN.resize(0);
	// Find min position of nodes which is infuenced by this particle
	Vector3i minx 	= Vector3i::Zero();
	Vector3i maxx 	= Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		minx(d) = (int) trunc(p0->X(d) - p0->PSize(d)-1.);
		maxx(d) = (int) ceil(p0->X(d) + p0->PSize(d)+1.);
	}

	// Find nodes within the influence range
	for (int i=minx(0); i<=maxx(0); ++i)
	for (int j=minx(1); j<=maxx(1); ++j)
	for (int k=minx(2); k<=maxx(2); ++k)
	{
		double n;
		Vector3d gn;
		Vector3d xn (i,j,k);
		NGN(p0->X, xn, Dx, p0->PSize, n, gn);
		// Find id of current node
		int id = i+j*Ncy+k*Ncz;
		if (n>0.)
		{
			p0->Lni.push_back(id);
			p0->LnN.push_back(n);
			p0->LnGN.push_back(gn);
		}
		// #pragma omp critical
		// {
		// 	Ln[id]->MPs.push_back(p0->ID);
		// }
	}
}

void MPM::ParticleToNode()
{
	// reset mass internal force velocity for nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		Ln[id]->Reset();
    }
    LAn.resize(0);
    // Update shape function and grad
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	CalNGN(Lp[p]);
    	Matrix3d vsp = -Lp[p]->Vol*Lp[p]->Stress;
		Vector3d fex = Lp[p]->M*Lp[p]->B + Lp[p]->Fh + Lp[p]->Fc;

		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 		n 	= Lp[p]->LnN[l];
			Vector3d 	gn 	= Lp[p]->LnGN[l];
			Vector3d 	df 	= n*fex + vsp*gn;
			// weigthed mass contribution
			double nm = n*Lp[p]->M;
			if (nm<0.)
			{
				cout << "mass problem" << endl;
				cout << "nm= " << nm << endl;
				cout << "p= " << p << endl;
				cout << "id= " << id << endl;
				cout << "Lp[p]->X= " << Lp[p]->X.transpose() << endl;
				cout << "Ln[id]->X= " << Ln[id]->X.transpose() << endl;
				abort();
			}
			#pragma omp atomic
			Ln[id]->M += nm;
			// #pragma omp atomic
			// Ln[id]->Mv(0) += nm*Lp[p]->V(0);
			// #pragma omp atomic
			// Ln[id]->Mv(1) += nm*Lp[p]->V(1);
			// #pragma omp atomic
			// Ln[id]->Mv(2) += nm*Lp[p]->V(2);
			// #pragma omp atomic
			// Ln[id]->F(0) += df(0);
			// #pragma omp atomic
			// Ln[id]->F(1) += df(1);
			// #pragma omp atomic
			// Ln[id]->F(2) += df(2);
			// #pragma omp atomic
			// Ln[id]->Mv(0) += df(0)*Dt;
			// #pragma omp atomic
			// Ln[id]->Mv(1) += df(1)*Dt;
			// #pragma omp atomic
			// Ln[id]->Mv(2) += df(2)*Dt;
			for (size_t d=0; d<D; ++d)
			{
				#pragma omp atomic
				Ln[id]->Mv(d) += nm*Lp[p]->V(d);
				#pragma omp atomic
				Ln[id]->F(d) += df(d);
				#pragma omp atomic
				Ln[id]->Mv(d) += df(d)*Dt;
				// smooth stress on node
				for (size_t c=0; c<D; ++c)
				{
					#pragma omp atomic
					Ln[id]->Stress(d,c) += nm*Lp[p]->Stress(d,c);
				}
			}
		}
    }
    UpdateLAn();
}

void MPM::CalNGN_MLS(MPM_PARTICLE* p0)
{
	// Reset shape function (N) and gradient of shape function (GN)
	p0->Lni.resize(0);
	p0->Lgi.resize(0);
	p0->LnN.resize(0);
	p0->LnGN.resize(0);
	// Find min position of nodes which is infuenced by this particle
	// for nodes
	Vector3i minn 	= Vector3i::Zero();
	Vector3i maxn 	= Vector3i::Zero();
	// for gauss points
	Vector3i ming 	= Vector3i::Zero();
	Vector3i maxg 	= Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		// minn(d) = ceil(p0->X(d) -1.);
		// maxn(d) = trunc(p0->X(d) +1.);
		minn(d) = ceil(p0->X(d) -1.5);
		maxn(d) = trunc(p0->X(d) +1.5);
		ming(d) = ceil(p0->X(d)-2.);
		maxg(d) = trunc(p0->X(d)+1.);
	}

	// Find nodes within the influence range
	for (int i=minn(0); i<=maxn(0); ++i)
	for (int j=minn(1); j<=maxn(1); ++j)
	for (int k=minn(2); k<=maxn(2); ++k)
	{
		double n;
		Vector3d gn;
		// Find id of current node
		int id = i+j*Ncy+k*Ncz;
		if (id>Ln.size())
		{
			cout << "id too large" << endl;
			cout << p0->X.transpose() << endl;
			cout << i << " " << j << " " << k << endl;
			cout << id << endl;
			cout << p0->ID << endl;
			cout << Lp.size() << endl;
		}
		NGN(p0->X, Ln[id]->X, Dx, p0->PSize, n, gn);
		if (n>0.)
		{
			p0->Lni.push_back(id);
			p0->LnN.push_back(n);
			p0->LnGN.push_back(gn);
			Ln[id]->Actived = true;
		}
		// store particle ID for MLS
		#pragma omp critical
		{
			Ln[id]->MPs.push_back(p0->ID);
		}
	}
	// Find nodes within the influence range
	// for (int i=ming(0); i<=maxg(0); ++i)
	// for (int j=ming(1); j<=maxg(1); ++j)
	// for (int k=ming(2); k<=maxg(2); ++k)
	// {
	// 	// Find id of current node
	// 	int id = i+j*Ncy+k*Ncz;
	// 	if (id>Ln.size())
	// 	{
	// 		cout << "id too large GP" << endl;
	// 		cout << p0->X.transpose() << endl;
	// 		cout << i << " " << j << " " << k << endl;
	// 		cout << id << endl;
	// 		cout << p0->ID << endl;
	// 	}
	// 	// store particle ID for MLS
	// 	#pragma omp critical
	// 	{
	// 		Lg[id]->MPs.push_back(p0->ID);
	// 		// if (Lg[id]->MPs.size()>50)
	// 		// {
	// 		// 	cout << "Lg[id]->MPs too large" << endl;
	// 		// 	cout << Lg[id]->MPs.size() << endl;
	// 		// 	abort();
	// 		// }
	// 	}
	// 	p0->Lgi.push_back(id);
	// }
}

void MPM::ParticleToNodeMLS()
{
	// reset mass internal force velocity for nodes
	// #pragma omp parallel for schedule(static) num_threads(1)
	// for (size_t n=0; n<LAn.size(); ++n)
	// {
	// 	size_t id = LAn[n];
	// 	Ln[id]->Reset();
 //    }
    for (size_t n=0; n<Ln.size(); ++n)
    {
    	Ln[n]->Reset();
    }
    LAn.resize(0);

    // Update shape function and grad
    // #pragma omp parallel for schedule(static) num_threads(1)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	// UpdateLn(Lp[p]);
    	CalNGN_MLS(Lp[p]);

    	Matrix3d vsp = -Lp[p]->Vol*Lp[p]->Stress;
		Vector3d fex = Lp[p]->M*Lp[p]->B + Lp[p]->Fh;

		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
		{
			// Grid id
			size_t id = Lp[p]->Lni[l];
			// weight
			double 		n 	= Lp[p]->LnN[l];
			Vector3d 	gn 	= Lp[p]->LnGN[l];
			Vector3d 	df 	= n*fex + vsp*gn;
			// weigthed mass contribution
			double nm = n*Lp[p]->M;

			if (nm<0.)
			{
				cout << "mass problem" << endl;
				cout << "nm= " << nm << endl;
				cout << "p= " << p << endl;
				cout << "id= " << id << endl;
				cout << "Lp[p]->X= " << Lp[p]->X.transpose() << endl;
				cout << "Ln[id]->X= " << Ln[id]->X.transpose() << endl;
				abort();
			}
			#pragma omp atomic
			Ln[id]->M += nm;
			#pragma omp atomic
			Ln[id]->F(0) += df(0);
			#pragma omp atomic
			Ln[id]->F(1) += df(1);
			#pragma omp atomic
			Ln[id]->F(2) += df(2);
		}
    }
    UpdateLAn();
    // cout << "start mls" << endl;
	// #pragma omp parallel for schedule(static) num_threads(1)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		size_t nmp = Ln[id]->MPs.size();
		bool wrong = false;
		if (nmp>800)
		{
			cout << "nmp: " << nmp << endl;
			wrong = true;
		}
		// position of MPs
		vector<Vector3d> Xp;
		// velocity of MPs
		MatrixXd Vp(nmp,3);
		for (size_t i=0; i<nmp; ++i)
		{
			// ID of MP
			size_t p = Ln[id]->MPs[i];
			Xp.push_back(Lp[p]->X);
			Vp.row(i) = Lp[p]->V.transpose();
			if (wrong)
			{
				cout << p << endl;
			}
		}
		// shape functions
		// VectorXd phi = MLS(Xp, Ln[id]->X, 1);
		VectorXd phi = MLS(Xp, Ln[id]->X, 0);
		// // node velocity
		Ln[id]->V = phi.transpose()*Vp;
		Ln[id]->Mv = Ln[id]->M*Ln[id]->V+Ln[id]->F;
	}
	// cout << "end mls" << endl;
}

void MPM::UpdateLAn()
{
	vector<vector <size_t>> lan;
	lan.resize(Nproc);
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		auto id = omp_get_thread_num();
		lan[id].insert( lan[id].end(), Lp[p]->Lni.begin(), Lp[p]->Lni.end() );
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Nproc; ++n)
	{
		sort( lan[n].begin(), lan[n].end() );
		lan[n].erase(unique(lan[n].begin(), lan[n].end()), lan[n].end());
	}
	for (size_t n=0; n<Nproc; ++n)
	{
		LAn.insert( LAn.end(), lan[n].begin(), lan[n].end() );
	}
	sort( LAn.begin(), LAn.end() );
	LAn.erase(unique(LAn.begin(), LAn.end()), LAn.end());
}

void MPM::CalVOnNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		if (Ln[id]->M==0.)	Ln[id]->V.setZero();
		else
		{
			Vector3d fdamp = Dc*Ln[id]->F.norm()*Ln[id]->Mv.normalized();
			// cout << "fdamp= " << fdamp.transpose() << endl;
			Ln[id]->F  -= fdamp;
			Ln[id]->Mv -= fdamp*Dt;
			// Apply slipping BC
			if (Ln[id]->BCTypes.size()>0)
			{
				for (size_t i=0; i<Ln[id]->BCTypes.size(); ++i)
				{
					if (Ln[id]->BCTypes[i]==1)			Ln[id]->NonSlippingBC();
					else if (Ln[id]->BCTypes[i]==2)		Ln[id]->SlippingBC(Ln[id]->Norms[i]);
					else if (Ln[id]->BCTypes[i]==3)		Ln[id]->FrictionBC(Dt, Ln[id]->Norms[i]);
				}
			}
			Ln[id]->V = Ln[id]->Mv/Ln[id]->M;
			Ln[id]->Stress /= Ln[id]->M;
			// if (Ln[id]->M<1.0e-6)
			// {
			// 	// cout << "very small node mass" << endl;
			// 	Ln[id]->V.setZero();
			// }
		}
	}
}
void MPM::SetNonSlippingBC(size_t n)
{
	Ln[n]->BCTypes.push_back(1);
}
void MPM::SetNonSlippingBC(size_t i, size_t j, size_t k)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(1);
}

void MPM::SetSlippingBC(size_t n, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(2);
	Ln[n]->Norms.push_back(norm);
}
void MPM::SetSlippingBC(size_t i, size_t j, size_t k, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	// cout << "n= " << n << endl;
	// cout << "Nnode= " << Nnode << endl;
	Ln[n]->BCTypes.push_back(2);
	// cout << "push_back 1 " << endl;
	Ln[n]->Norms.push_back(norm);
	// cout << "push_back 2 " << endl;
}

void MPM::SetFrictionBC(size_t n, double mu, Vector3d& norm)
{
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;

}
void MPM::SetFrictionBC(size_t i, size_t j, size_t k, double mu, Vector3d& norm)
{
	int n = i+j*Ncy+k*Ncz;
	Ln[n]->BCTypes.push_back(3);
	Ln[n]->Norms.push_back(norm);
	Ln[n]->Mu = mu;
}

void MPM::NodeToParticle()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Reset position increasement of this particle
		Lp[p]->DeltaX 	= Vector3d::Zero();
		Lp[p]->StressSmooth.setZero();
 		if (!Lp[p]->FixV)
		{
			for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
			{
				size_t 	id = Lp[p]->Lni[l];
				double 	n  = Lp[p]->LnN[l];
				// Update velocity of this particle
				Vector3d an = Ln[id]->F/Ln[id]->M;
				Lp[p]->V += n*an*Dt;
				Lp[p]->X += n*Ln[id]->V*Dt;
				Lp[p]->StressSmooth += n*Ln[id]->Stress;
				// Lp[p]->X += n*Ln[id]->Mv/Ln[id]->M*Dt;
			}
		}
		else
		{
			Lp[p]->V = Lp[p]->Vf;
			Lp[p]->X += Lp[p]->V*Dt;
		}
		// if (Lp[p]->V.norm()>0.1)
		// {
		// 	cout << "p= " << p << endl;
		// 	cout << "Lp[p]->V: " << Lp[p]->V.transpose() << endl;
		// 	cout << "Lp[p]->X: " << Lp[p]->X.transpose() << endl;
		// 	abort();
		// }
		// Velocity gradient tensor
		CalVGradLocal(p);
		// Update deformation tensor
		Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->F;
		// Update particle length
		if (Lp[p]->Type==0)	CalPSizeR(p);
		// CalPSizeCP(p);
		// Update volume of particles
		Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
		// Update stress
		Matrix3d w = 0.5*Dt*((Lp[p]->L - Lp[p]->L.transpose()));
		Lp[p]->Stress += w*Lp[p]->Stress-Lp[p]->Stress*w.transpose();
		if (Lp[p]->Type==0)			Lp[p]->Elastic(de);
		else if (Lp[p]->Type==1)
		{
			Lp[p]->EOSMonaghan(Cs);
			Lp[p]->Newtonian(de);
		}
		else if (Lp[p]->Type==2)	Lp[p]->MohrCoulomb(de);
		else if (Lp[p]->Type==3)	Lp[p]->DruckerPrager(de);
		else if (Lp[p]->Type==5)
		{
			// Lp[p]->EOSMonaghan(Cs);
			Lp[p]->EOSMorris(Cs);
			Lp[p]->Granular(de);
		}
		// Reset hydro force and contact force
		Lp[p]->Fh.setZero();
		Lp[p]->Fc.setZero();
	}
	// int npf = 0;
	// for (size_t p=0; p<Lp.size(); ++p)
	// {
	// 	if (Lp[p]->Type==1)
	// 	{
	// 		npf++;
	// 	}
	// }
	// if (npf!=16200)
	// {
	// 	cout << "npf= " << npf << endl;
	// 	abort();
	// }
}

// void MPM::SmoothStress()
// {}

// void MPM::CalVOnNodeDoubleMapping()
// {
// 	#pragma omp parallel for schedule(static) num_threads(1)
// 	for (size_t c=0; c<LAn.size(); ++c)
// 	{
// 		int i = LAn[c][0];
// 		int j = LAn[c][1];
// 		int k = LAn[c][2];
// 		V[i][j][k] = Vector3d::Zero();
// 	}
// 	// Double map velocity from particles to nodes to aviod small mass problem
// 	#pragma omp parallel for schedule(static) num_threads(1)
// 	for (size_t p=0; p<Lp.size(); ++p)
// 	{
// 		for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
// 		{
// 			Vector3i ind = Lp[p]->Lni[l];
// 			double n = Lp[p]->LnN[l];
// 			// Update velocity of this node
// 			#pragma omp critical
// 			{
// 				V[ind(0)][ind(1)][ind(2)] += n*Lp[p]->M*Lp[p]->V/M[ind(0)][ind(1)][ind(2)];		
// 			}
// 		}
// 	}
// }

void MPM::CalVGradLocal(int p)
{
	Lp[p]->L = Matrix3d::Zero();
	for (size_t l=0; l<Lp[p]->Lni.size(); ++l)
	{
		size_t	 	id = Lp[p]->Lni[l];
		Vector3d 	gn 	= Lp[p]->LnGN[l];
		// Calculate velocity gradient tensor
		Lp[p]->L += gn*Ln[id]->V.transpose();
	}
}

void MPM::CalPSizeCP(int p)
{
	Lp[p]->PSize(0) =  Lp[p]->PSize0(0)*Lp[p]->F(0,0);
	Lp[p]->PSize(1) =  Lp[p]->PSize0(1)*Lp[p]->F(1,1);
	Lp[p]->PSize(2) =  Lp[p]->PSize0(2)*Lp[p]->F(2,2);
}

// Based on "iGIMP: An implicit generalised interpolation material point method for large deformations"
void MPM::CalPSizeR(int p)
{
	Lp[p]->PSize(0) =  Lp[p]->PSize0(0)*sqrt(Lp[p]->F(0,0)*Lp[p]->F(0,0) + Lp[p]->F(1,0)*Lp[p]->F(1,0) + Lp[p]->F(2,0)*Lp[p]->F(2,0));
	Lp[p]->PSize(1) =  Lp[p]->PSize0(1)*sqrt(Lp[p]->F(0,1)*Lp[p]->F(0,1) + Lp[p]->F(1,1)*Lp[p]->F(1,1) + Lp[p]->F(2,1)*Lp[p]->F(2,1));
	Lp[p]->PSize(2) =  Lp[p]->PSize0(2)*sqrt(Lp[p]->F(0,2)*Lp[p]->F(0,2) + Lp[p]->F(1,2)*Lp[p]->F(1,2) + Lp[p]->F(2,2)*Lp[p]->F(2,2));
}

void MPM::CalStressOnParticleElastic()
{
	// Update stresses on particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Velocity gradient tensor
		CalVGradLocal(p);
		// Update deformation tensor
		Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->F;
		// Update particle length
		CalPSizeR(p);
		// CalPSizeCP(p);
		// Update volume of particles
		Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
		// Update stress
		Matrix3d w = 0.5*Dt*((Lp[p]->L - Lp[p]->L.transpose()));
		Lp[p]->Stress += w*Lp[p]->Stress-Lp[p]->Stress*w.transpose();
		Lp[p]->Elastic(de);
	}
}

void MPM::CalStressOnParticleMohrCoulomb()
{
	// Update stresses on particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Velocity gradient tensor
		CalVGradLocal(p);
		// Update deformation tensor
		Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L)*Lp[p]->F;
		// Update particle length
		// CalPSizeR(p);
		CalPSizeCP(p);
		// Update volume of particles
		Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
		// Update stress
		Lp[p]->MohrCoulomb(de);
		// if (Lp[p]->Stress(2,2)<-0.02)
		// {
		// 	cout << "Lp[p]->Stress(2,2)<-0.02" << endl;
		// 	cout << "p= " << p << endl;
		// 	abort();
		// }
	}
}

void MPM::CalStressOnParticleNewtonian()
{
	// cout << "start CalStressOnParticleNewtonian " << endl;
	// Update stresses on particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Velocity gradient tensor
		CalVGradLocal(p);
		// Update deformation tensor
		Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L*Dt)*Lp[p]->F;
		// Update particle length
		// CalPSizeCP(p);
		// Update volume of particles
		Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*Dt*(Lp[p]->L + Lp[p]->L.transpose());
		// Update EOS
		// cout << "start EOSMorris " << endl;
		// Lp[p]->EOSMorris(Cs);
		Lp[p]->EOSMonaghan(Cs);
		// Update stress
		// cout << "start Newtonian " << endl;
		Lp[p]->Newtonian(de);
	// cout << "finish CalStressOnParticleNewtonian " << endl;

	}
}

void MPM::SolveMUSL(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		bool show = false;
		if (t%100==0)	show = true;
		if (show) 	cout << "Time Step = " << t << endl;
		if (t%ts == 0)
		// if (t> 80600)
		{
			WriteFileH5(t);
		}
		auto t_start = std::chrono::system_clock::now();
		if (MLSv)	ParticleToNodeMLS();
		else 		ParticleToNode();
		auto t_end = std::chrono::system_clock::now();
		if (show)	cout << "ParticleToNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		CalVOnNode();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalVOnNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		// if (t==2)
		// {
		// 	cout << Lp[4]->P << endl;
		// }
		// cout << "t= " << t << endl;
		// cout << Lp[4]->P << endl;
		NodeToParticle();
		// cout << "t= " << t << endl;
		// cout << Lp[4]->P << endl;

		t_end = std::chrono::system_clock::now();
		if (show)	cout << "NodeToParticle= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		// if 		(CMType==0) 	CalStressOnParticleElastic();
		// else if (CMType==1) 	CalStressOnParticleMohrCoulomb();
		// else if (CMType==2) 	CalStressOnParticleNewtonian();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalStressOnParticleMohrCoulomb= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (show) 	cout << "===========================" << endl;
	}
}

// void MPM::SolveUSF(int tt, int ts)
// {
// 	for (int t=0; t<tt; ++t)
// 	{
// 		if (t%ts == 0)
// 		{
// 			cout << "Time Step = " << t << endl;
// 			WriteFileH5(t);
// 		}

// 		ParticleToNode();
// 		CalVOnNode();

// 		for (size_t i=0; i<LFn.size(); ++i)
// 		{
// 			V[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 		}

// 		if 		(CMType==0) 	CalStressOnParticleElastic();
// 		else if (CMType==1) 	CalStressOnParticleMohrCoulomb();
// 		CalFOnNode(false);

// 		for (size_t i=0; i<LFn.size(); ++i)
// 		{
// 			F[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 			Mv[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 		}

// 		NodeToParticle();
// 	}
// }

// void MPM::SolveUSA(int tt, int ts)
// {
// 	for (int t=0; t<tt; ++t)
// 	{
// 		if (t%ts == 0)
// 		{
// 			cout << "Time Step = " << t << endl;
// 			WriteFileH5(t);
// 		}

// 		ParticleToNode();
// 		CalVOnNode();

// 		for (size_t i=0; i<LFn.size(); ++i)
// 		{
// 			V[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 		}

// 		#pragma omp parallel for schedule(static) num_threads(1)
// 		for (size_t p=0; p<Lp.size(); ++p)
// 		{
// 			CalVGradLocal(p);
// 			// Update strain
// 			Matrix3d de = 0.25*(Lp[p]->L + Lp[p]->L.transpose());
// 			// Update stress
// 			Lp[p]->Stress += 2.*Lp[p]->Mu*de + Lp[p]->La*de.trace()*Matrix3d::Identity();
// 		}

// 		CalFOnNode(false);

// 		for (size_t i=0; i<LFn.size(); ++i)
// 		{
// 			F[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 			Mv[LFn[i][0]][LFn[i][1]][LFn[i][2]](LFn[i][3]) = 0.;
// 		}

// 		NodeToParticle();
// 		CalVOnNodeDoubleMapping();

// 		#pragma omp parallel for schedule(static) num_threads(1)
// 		for (size_t p=0; p<Lp.size(); ++p)
// 		{	
// 			Matrix3d L0 = Lp[p]->L;
// 			CalVGradLocal(p);
// 			// Update strain
// 			Matrix3d de = 0.25*(Lp[p]->L + Lp[p]->L.transpose());
// 			// Update stress
// 			Lp[p]->Stress += 2.*Lp[p]->Mu*de + Lp[p]->La*de.trace()*Matrix3d::Identity();
// 			// Update deformation tensor
// 			Lp[p]->F = (Matrix3d::Identity() + 0.5*(L0+Lp[p]->L))*Lp[p]->F;
// 			// Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L)*Lp[p]->F;
// 			// Update particle length
// 			CalPSizeR(p);
// 			// Update volume of particles
// 			Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
// 		}
// 	}
// }

void MPM::AddNode(size_t level, Vector3d& x)
{
    Ln.push_back(new MPM_NODE(level,x));
    Ln[Ln.size()-1]->ID = Ln.size()-1;
}

// void MPM::AddParticle(int tag, Vector3d& x, double m, double young, double poisson)
// {
//     Lp.push_back(new MPM_PARTICLE(tag,x,m, young, poisson));
//     Lp[Lp.size()-1]->ID = Lp.size()-1;
// }

void MPM::AddParticle(int tag, Vector3d& x, double m)
{
    Lp.push_back(new MPM_PARTICLE(tag,x,m));
    Lp[Lp.size()-1]->ID = Lp.size()-1;
}

void MPM::DeleteParticles()
{
	vector <MPM_PARTICLE*>	Lpt;
	Lpt.resize(0);

	for (size_t p=0; p<Lp.size(); ++p)
	{
		if (!Lp[p]->Removed)	Lpt.push_back(Lp[p]);
	}
	Lp = Lpt;

	for (size_t p=0; p<Lp.size(); ++p)
	{
		Lp[p]->ID = p;
	}
}

// void MPM::AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson)
// {
// 	Vector3i maxx = Vector3i::Zero();
// 	for (size_t d=0; d<D; ++d)
// 	{
// 		maxx(d) = (int) (l(d)/ratio)-1;
// 	}

// 	for (int k=0; k<=maxx(2); ++k)
//     for (int j=0; j<=maxx(1); ++j)
//     for (int i=0; i<=maxx(0); ++i)
//     {
//     	int type = -1;

//     	if (D==1)
//     	{
//     		if (i==0 || i==maxx(0))	type = -2;
//     	}
//     	else if (D==2)
//     	{
//     		if (i==0 || i==maxx(0) || j==0 || j==maxx(1))	type = -2;		
//     	}
//     	else if (D==3)
//     	{
//     		if (i==0 || i==maxx(0) || j==0 || j==maxx(1) || k==0 || k==maxx(2))		type = -2;
//     	}

//     	Vector3d x = Vector3d::Zero();
//     					x(0) = ratio*(i + 0.5)+x0(0);
//     	if (D>1)		x(1) = ratio*(j + 0.5)+x0(1);
//     	if (D>2)		x(2) = ratio*(k + 0.5)+x0(2);

//     	AddParticle(type, x, m, young, poisson);
//     }

//     Lbp.resize(0);
//     for (size_t p=0; p<Lp.size(); ++p)
//     {
//     	Lp[p]->Vol0 = 1.;
//     	for (size_t d=0; d<D; ++d)
//     	{
//     		Lp[p]->Vol0 *= ratio;
//     		if (Ntype==3)	Lp[p]->PSize0(d) = 0.5*ratio;
//     		else 			Lp[p]->PSize0(d) = 0.;
//     		Lp[p]->PSize(d) = Lp[p]->PSize0(d);
//     	}
//     	Lp[p]->Vol 	= Lp[p]->Vol0;

//     	if (Lp[p]->Type==-2)
//     	{
//     		Lbp.push_back(Lp[p]);
//     	}
//     }
// }

void MPM::AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m)
{
	Vector3i maxx = Vector3i::Zero();
	for (size_t d=0; d<D; ++d)
	{
		maxx(d) = (int) (l(d)/ratio)-1;
	}

	for (int k=0; k<=maxx(2); ++k)
    for (int j=0; j<=maxx(1); ++j)
    for (int i=0; i<=maxx(0); ++i)
    {
    	Vector3d x = Vector3d::Zero();
    					x(0) = ratio*(i + 0.5)+x0(0);
    	if (D>1)		x(1) = ratio*(j + 0.5)+x0(1);
    	if (D>2)		x(2) = ratio*(k + 0.5)+x0(2);

    	AddParticle(tag, x, m);
    }

    for (size_t p=0; p<Lp.size(); ++p)
    {
    	Lp[p]->Vol0 = 1.;
    	for (size_t d=0; d<D; ++d)
    	{
    		Lp[p]->Vol0 *= ratio;
    		if (Ntype==3)	Lp[p]->PSize0(d) = 0.5*ratio;
    		else 			Lp[p]->PSize0(d) = 0.;
    		Lp[p]->PSize(d) = Lp[p]->PSize0(d);
    	}
    	Lp[p]->Vol 	= Lp[p]->Vol0;
    	// Assume a radius of MPM particle for DEMPM 
    	Lp[p]->R 	= pow(3.*Lp[p]->Vol0/(4.*M_PI), 1./3.);
    }
}

inline void MPM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "MPM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};			//create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.
	hsize_t	dims_tensor[1] = {6*Lp.size()};

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);
	int rank_tensor = sizeof(dims_tensor) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);
	DataSpace	*space_tensor = new DataSpace(rank_tensor, dims_tensor);

	double* tag_h5 	= new double[  Lp.size()];
	double* m_h5 	= new double[  Lp.size()];
	double* you_h5 	= new double[  Lp.size()];
	double* poi_h5 	= new double[  Lp.size()];
	double* pos_h5 	= new double[3*Lp.size()];
	double* vel_h5 	= new double[3*Lp.size()];
	double* s_h5 	= new double[6*Lp.size()];

	double* szz_h5 	= new double[  Lp.size()];

	for (size_t i=0; i<Lp.size(); ++i)
	{
        tag_h5[  i  ] 	= Lp[i]->Tag;
        m_h5  [  i  ] 	= Lp[i]->M;
        you_h5[  i  ] 	= Lp[i]->Young;
        poi_h5[  i  ] 	= Lp[i]->Poisson;
		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);
		s_h5  [6*i  ] 	= Lp[i]->StressSmooth(0,0);
		s_h5  [6*i+1] 	= Lp[i]->StressSmooth(0,1);
		s_h5  [6*i+2] 	= Lp[i]->StressSmooth(0,2);
		s_h5  [6*i+3] 	= Lp[i]->StressSmooth(1,1);
		s_h5  [6*i+4] 	= Lp[i]->StressSmooth(1,2);
		s_h5  [6*i+5] 	= Lp[i]->StressSmooth(2,2);

		// szz_h5[i] 		= Lp[i]->StressSmooth(1,1);
		szz_h5[i] 		= Lp[i]->P;
		// szz_h5[i] 		= Lp[i]->Stress(1,1);
	}

	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_m		= new DataSet(file.createDataSet("Mass", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_you	= new DataSet(file.createDataSet("Young", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_poi	= new DataSet(file.createDataSet("Poisson", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_s		= new DataSet(file.createDataSet("Stress", PredType::NATIVE_DOUBLE, *space_tensor));

    DataSet	*dataset_szz		= new DataSet(file.createDataSet("SZZ", PredType::NATIVE_DOUBLE, *space_scalar));

	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
	dataset_m->write(m_h5, PredType::NATIVE_DOUBLE);
	dataset_you->write(you_h5, PredType::NATIVE_DOUBLE);
	dataset_poi->write(poi_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
	dataset_s->write(s_h5, PredType::NATIVE_DOUBLE);

	dataset_szz->write(szz_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;
	delete space_tensor;
	delete dataset_tag;
	delete dataset_m;
	delete dataset_you;
	delete dataset_poi;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_s;

	delete dataset_szz;

	delete tag_h5;
	delete m_h5;
	delete you_h5;
	delete poi_h5;
	delete pos_h5;
	delete vel_h5;
	delete s_h5;

	delete szz_h5;

	file.close();

	string file_name_xmf = "MPM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Lp.size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Lp.size() << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"SZZ\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/SZZ \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Stress\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}
