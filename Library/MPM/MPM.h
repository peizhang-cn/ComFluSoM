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

#include "./HEADER.h"
#include <SHAPE.h>
#include <MPM_PARTICLE.h>

class MPM
{
public:
	MPM();
	~MPM();
	MPM(int ntype, int nx, int ny, int nz, Vector3d dx);
	void Init();
	void UpdateLn(MPM_PARTICLE* p0);
	void CalNGN(MPM_PARTICLE* p0);
	void ParticleToNode();
	void CalFOnNode(bool firstStep);
	void NodeToParticle();
	void CalStressOnParticle();
	void CalVGradLocal(int p);
	void CalPSizeCP(int p);
	void CalPSizeR(int p);
	void CalVOnNode();
	void CalVOnNodeDoubleMapping();
	void ApplyFrictionBC(int i, int j, int k, Vector3d n, double muf);
	void SolveMUSL(int tt, int ts);
	void SolveUSF(int tt, int ts);
	void SolveUSA(int tt, int ts);
	void AddParticle(int tag, Vector3d& x, double m, double young, double poisson);
	void DeleteParticles();
	void AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson);
	void WriteFileH5(int n);

	double 		(*N)(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp);
	Vector3d 	(*GN)(Vector3d x, Vector3d xc, Vector3d l, Vector3d lp);

	double*** 						M;														// Mass of nodes
	Vector3d***						V;														// Velocity of nodes
	Vector3d***						Mv;														// Momentum of nodes
	Vector3d***						F;														// Total force of nodes

	vector <Vector3i>				LAn;													// List of actived nodes
	vector <Vector3i>				LFn;													// List of nodes for applying fixed boundary condition

	vector <MPM_PARTICLE*>			Lp;														// List of all MPM particles
	vector <MPM_PARTICLE*>			Lbp;													// List of boundary MPM particles

    int 							Nx;														// Domain size
    int 							Ny;
    int 							Nz;

    int 							Nproc;

    int 							D;														// Dimension	
    int 							Ntype;													// Type of shape function 0 for Linear, 1 for Quadratic and 2 for Cubic 3 for GIMP
    double 							Nrange;													// Influence range of shape function

    Vector3d						Dx;														// Space step
};

MPM::MPM(int ntype, int nx, int ny, int nz, Vector3d dx)
{
	Nproc	= 1;

	Ntype 	= ntype;
	Nx 		= nx;
	Ny 		= ny;
	Nz 		= nz;
	Dx 		= dx;
	D 		= 3;

	if (Nz==0)
	{
		D = 2;
		if (Ny==0)	D = 1;
	}
	// Linear
	if 		(Ntype == 0)
	{
		if (D==1)
		{
			N  		=& ShapeL1D;
			GN 		=& GradShapeL1D;			
		}
		else if (D==2)
		{
			N  		=& ShapeL2D;
			GN 		=& GradShapeL2D;
		}
		else if (D==3)
		{
			N  		=& ShapeL3D;
			GN 		=& GradShapeL3D;
		}
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.;
		cout << "Using Linear shape function." << endl;
	}
	// Quadratic
	else if (Ntype == 1)
	{
		if (D==1)
		{
			N  		=& ShapeQ1D;
			GN 		=& GradShapeQ1D;			
		}
		else if (D==2)
		{
			N  		=& ShapeQ2D;
			GN 		=& GradShapeL2D;
		}
		else if (D==3)
		{
			N  		=& ShapeQ3D;
			GN 		=& GradShapeQ3D;
		}
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 1.5;
		cout << "Using Quadratic shape function." << endl;
	}
	// Cubic
	else if (Ntype == 2)
	{
		if (D==1)
		{
			N  		=& ShapeC1D;
			GN 		=& GradShapeC1D;			
		}
		else if (D==2)
		{
			N  		=& ShapeC2D;
			GN 		=& GradShapeC2D;
		}
		else if (D==3)
		{
			N  		=& ShapeC3D;
			GN 		=& GradShapeC3D;
		}
		else
		{
			cout << "Dimension is higher than 3." << endl;
			abort();
		}
		Nrange 	= 2.;
		cout << "Using Cubic shape function." << endl;
	}
	// GIMP
	else if (Ntype == 3)
	{
		if (D==1)
		{
			N  		=& ShapeGIMP1D;
			GN 		=& GradShapeGIMP1D;			
		}
		else if (D==2)
		{
			N  		=& ShapeGIMP2D;
			GN 		=& GradShapeGIMP2D;
		}
		else if (D==3)
		{
			N  		=& ShapeGIMP3D;
			GN 		=& GradShapeGIMP3D;
		}
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
	LFn.resize(0);

	M	= new double** 		[Nx+1];
	V	= new Vector3d** 	[Nx+1];
	Mv	= new Vector3d** 	[Nx+1];
	F	= new Vector3d** 	[Nx+1];

	for (int i=0; i<=Nx; ++i)
	{
		M [i]	= new double* 	[Ny+1];
		V [i]	= new Vector3d* [Ny+1];
		Mv[i]	= new Vector3d* [Ny+1];
		F [i]	= new Vector3d* [Ny+1];

		for (int j=0; j<=Ny; ++j)
		{
			M [i][j]	= new double 	[Nz+1];
			V [i][j]	= new Vector3d 	[Nz+1];
			Mv[i][j]	= new Vector3d 	[Nz+1];
			F [i][j]	= new Vector3d 	[Nz+1];

			for (int k=0; k<=Nz; ++k)
			{
				M [i][j][k]		= 0.;
				V [i][j][k]		= Vector3d::Zero();
				Mv[i][j][k]		= Vector3d::Zero();
				F [i][j][k]		= Vector3d::Zero();
			}		
		}
	}
	cout << "================ Finish init. ================" << endl;
}

void MPM::UpdateLn(MPM_PARTICLE* p0)
{
	p0->Ln.resize(0);
	// Find min position of nodes which is infuenced by this particle
	Vector3i minx 	= Vector3i::Zero();
	Vector3i maxx 	= Vector3i::Zero();
	for (int d=0; d<D; ++d)
	{
		minx(d) = (int) trunc(p0->X(d) - Nrange - p0->PSize(d) + 1.);
		maxx(d) = (int) ceil(p0->X(d) + Nrange + p0->PSize(d) - 1.);
	}
	// Find nodes within the influence range
	for (int i=minx(0); i<=maxx(0); ++i)
	for (int j=minx(1); j<=maxx(1); ++j)
	for (int k=minx(2); k<=maxx(2); ++k)
	{
		// Find position index of current node
		Vector3i ind (i,j,k);
		p0->Ln.push_back(ind);
	}
}

void MPM::CalNGN(MPM_PARTICLE* p0)
{
	// Reset shape function (N) and gradient of shape function (GN)
	p0->LnN.resize(0);
	p0->LnGN.resize(0);	
	for (size_t l=0; l<p0->Ln.size(); ++l)
	{
		// Calculate shape function (N) and gradient of shape function (GN)
		Vector3d xc;
		xc << p0->Ln[l](0), p0->Ln[l](1), p0->Ln[l](2);
		double n = N(p0->X, xc, Dx, p0->PSize);
		Vector3d gn = GN(p0->X, xc, Dx, p0->PSize);
		p0->LnN.push_back(n);
		p0->LnGN.push_back(gn);
	}
}

void MPM::ParticleToNode()
{
	// reset mass internal force velocity for nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<LAn.size(); ++c)
	{
		int i = LAn[c](0);
		int j = LAn[c](1);
		int k = LAn[c](2);

    	M [i][j][k] = 0.;
    	F [i][j][k] = Vector3d::Zero();
    	Mv[i][j][k] = Vector3d::Zero();
    }
    LAn.resize(0);

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	UpdateLn(Lp[p]);
    	CalNGN(Lp[p]);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		for (size_t l=0; l<Lp[p]->Ln.size(); ++l)
		{
			Vector3i ind = Lp[p]->Ln[l];
			double n = Lp[p]->LnN[l];
			// Update mass and momentum of this node
			#pragma omp critical
			{
				M [ind(0)][ind(1)][ind(2)] += n*Lp[p]->M;
				Mv[ind(0)][ind(1)][ind(2)] += n*Lp[p]->M*Lp[p]->V;	
				LAn.push_back(ind);
			}
		}
	}
}

void MPM::CalFOnNode(bool firstStep)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		for (size_t l=0; l<Lp[p]->Ln.size(); ++l)
		{
			Vector3i 	ind = Lp[p]->Ln[l];
			double 		n   = Lp[p]->LnN[l];
			Vector3d 	gn 	= Lp[p]->LnGN[l];
			// Update force and momentum of this node
			Matrix3d 	vsp = -Lp[p]->Vol*Lp[p]->S;
			Vector3d 	df  = n*(Lp[p]->M*Lp[p]->B + Lp[p]->Fh) + vsp.transpose()*gn;
			if (firstStep)	df *= 0.5;

			#pragma omp critical
			{
				F [ind(0)][ind(1)][ind(2)] += df;
				Mv[ind(0)][ind(1)][ind(2)] += df;
			}
		}	
	}
}

// void MPM::ApplyFrictionBC(int i, int j, int k, Vector3d n, double muf)
// {
// 	double fnNorm = n.dot(F[i][j][k]);			// normal force

// 	if (fnNorm>0.)
// 	{
// 		Vector3d ft = F[i][j][k]-fnNorm*n;		// tangential force

// 		Vector3d mvn = n.dot(Mv[i][j][k])*n;	// normal momentum
// 		Vector3d mvt = Mv[i][j][k]-mvn;			// tangential momentum
// 		// F[i][j][k] = fnNorm*n;
// 		// Mv[i][j][k] = fnNorm*n;
		
// 		// Vector3d mvt0 = Mv[i][j][k] - mvn*n;
// 		// Vector3d vt = (V[i][j][k]-n.dot(V[i][j][k])).normalized()
// 		Vector3d fft = -(V[i][j][k]-n.dot(V[i][j][k])*n).normalized()*muf*fnNorm;		// tangential friction force
		
// 		Vector3d mvt1 = mvt + fft;

// 		F[i][j][k] 	= Vector3d::Zero();
// 		Mv[i][j][k] = Vector3d::Zero();

// 		if (mvt.dot(mvt1)>0.)
// 		{
// 			F[i][j][k]  = ft+fft/*+fnNorm*n*/;
// 			Mv[i][j][k] = mvt1/*+mvn*/;
// 		}
// 		// else
// 		// {
// 		// 	F[i][j][k]  = fnNorm*n;
// 		// 	Mv[i][j][k] = mvn;			
// 		// }
// 		// else
// 		// {
// 		// 	F[i][j][k] 	= Vector3d::Zero();
// 		// 	Mv[i][j][k] = Vector3d::Zero();
// 		// }
// 	}
// }

void MPM::NodeToParticle()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Reset position increasement of this particle
		Lp[p]->DeltaX 	= Vector3d::Zero();

		if (!Lp[p]->FixV)
		{
			for (size_t l=0; l<Lp[p]->Ln.size(); ++l)
			{
				Vector3i 	ind = Lp[p]->Ln[l];
				double 		n   = Lp[p]->LnN[l];
				// Update velocity of this particle
				Vector3d deltaV = F[ind(0)][ind(1)][ind(2)]/M[ind(0)][ind(1)][ind(2)];
				Lp[p]->V += n*deltaV;
				// Update the position increasement of this particle
				Vector3d deltaX = Mv[ind(0)][ind(1)][ind(2)]/M[ind(0)][ind(1)][ind(2)];
				Lp[p]->X += n*deltaX;
			}	
		}
		else (Lp[p]->V = Lp[p]->Vf);
	}		
}

void MPM::CalVOnNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<LAn.size(); ++c)
	{
		int i = LAn[c](0);
		int j = LAn[c](1);
		int k = LAn[c](2);
		if (M[i][j][k]==0.)		V[i][j][k] = Vector3d::Zero();
		else					V[i][j][k] = Mv[i][j][k]/M[i][j][k];
	}
}

void MPM::CalVOnNodeDoubleMapping()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<LAn.size(); ++c)
	{
		int i = LAn[c](0);
		int j = LAn[c](1);
		int k = LAn[c](2);
		V[i][j][k] = Vector3d::Zero();
	}
	// Double map velocity from particles to nodes to aviod small mass problem
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		for (size_t l=0; l<Lp[p]->Ln.size(); ++l)
		{
			Vector3i ind = Lp[p]->Ln[l];
			double n = Lp[p]->LnN[l];
			// Update velocity of this node
			#pragma omp critical
			{
				V[ind(0)][ind(1)][ind(2)] += n*Lp[p]->M*Lp[p]->V/M[ind(0)][ind(1)][ind(2)];		
			}
		}
	}
}

void MPM::CalVGradLocal(int p)
{
	Lp[p]->L = Matrix3d::Zero();
	for (size_t l=0; l<Lp[p]->Ln.size(); ++l)
	{
		Vector3i 	ind = Lp[p]->Ln[l];
		Vector3d 	gn 	= Lp[p]->LnGN[l];
		// Calculate velocity gradient tensor
		Lp[p]->L += gn*V[ind(0)][ind(1)][ind(2)].transpose();
	}
}

void MPM::CalPSizeCP(int p)
{
	Lp[p]->PSize(0) =  Lp[p]->PSize0(0)*Lp[p]->F(0,0);
	Lp[p]->PSize(1) =  Lp[p]->PSize0(1)*Lp[p]->F(1,1);
	Lp[p]->PSize(2) =  Lp[p]->PSize0(2)*Lp[p]->F(2,2);
}

void MPM::CalPSizeR(int p)
{
	Lp[p]->PSize(0) =  Lp[p]->PSize0(0)*sqrt(Lp[p]->F(0,0)*Lp[p]->F(0,0) + Lp[p]->F(1,0)*Lp[p]->F(1,0) + Lp[p]->F(2,0)*Lp[p]->F(2,0));
	Lp[p]->PSize(1) =  Lp[p]->PSize0(1)*sqrt(Lp[p]->F(0,1)*Lp[p]->F(0,1) + Lp[p]->F(1,1)*Lp[p]->F(1,1) + Lp[p]->F(2,1)*Lp[p]->F(2,1));
	Lp[p]->PSize(2) =  Lp[p]->PSize0(2)*sqrt(Lp[p]->F(0,2)*Lp[p]->F(0,2) + Lp[p]->F(1,2)*Lp[p]->F(1,2) + Lp[p]->F(2,2)*Lp[p]->F(2,2));
}

void MPM::CalStressOnParticle()
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
		CalPSizeR(p);
		// Update volume of particles
		Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*(Lp[p]->L + Lp[p]->L.transpose());
		// Update stress
		Lp[p]->S += 2.*Lp[p]->Mu*de + Lp[p]->La*de.trace()*Matrix3d::Identity();
		// Lp[p]->MCModel();
	}
}

void MPM::SolveMUSL(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			WriteFileH5(t);
		}

		ParticleToNode();

		if (t==0)	CalFOnNode(true);
		else 		CalFOnNode(false);

		for (size_t i=0; i<LFn.size(); ++i)
		{
			F[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
			Mv[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
		}

		NodeToParticle();
		CalVOnNode();
		CalStressOnParticle();
	}
}

void MPM::SolveUSF(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			WriteFileH5(t);
		}

		ParticleToNode();
		CalVOnNode();

		for (size_t i=0; i<LFn.size(); ++i)
		{
			V[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
		}

		CalStressOnParticle();
		CalFOnNode(false);

		for (size_t i=0; i<LFn.size(); ++i)
		{
			F[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
			Mv[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
		}

		NodeToParticle();
	}
}

void MPM::SolveUSA(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			WriteFileH5(t);
		}

		ParticleToNode();
		CalVOnNode();

		for (size_t i=0; i<LFn.size(); ++i)
		{
			V[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
		}

		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t p=0; p<Lp.size(); ++p)
		{
			CalVGradLocal(p);
			// Update strain
			Matrix3d de = 0.25*(Lp[p]->L + Lp[p]->L.transpose());
			// Update stress
			Lp[p]->S += 2.*Lp[p]->Mu*de + Lp[p]->La*de.trace()*Matrix3d::Identity();
		}

		CalFOnNode(false);

		for (size_t i=0; i<LFn.size(); ++i)
		{
			F[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
			Mv[LFn[i][0]][LFn[i][1]][LFn[i][2]] = Vector3d::Zero();
		}

		NodeToParticle();
		CalVOnNodeDoubleMapping();

		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t p=0; p<Lp.size(); ++p)
		{	
			Matrix3d L0 = Lp[p]->L;
			CalVGradLocal(p);
			// Update strain
			Matrix3d de = 0.25*(Lp[p]->L + Lp[p]->L.transpose());
			// Update stress
			Lp[p]->S += 2.*Lp[p]->Mu*de + Lp[p]->La*de.trace()*Matrix3d::Identity();
			// Update deformation tensor
			Lp[p]->F = (Matrix3d::Identity() + 0.5*(L0+Lp[p]->L))*Lp[p]->F;
			// Lp[p]->F = (Matrix3d::Identity() + Lp[p]->L)*Lp[p]->F;
			// Update particle length
			CalPSizeR(p);
			// Update volume of particles
			Lp[p]->Vol 	= Lp[p]->F.determinant()*Lp[p]->Vol0;
		}
	}
}

void MPM::AddParticle(int tag, Vector3d& x, double m, double young, double poisson)
{
    Lp.push_back(new MPM_PARTICLE(tag,x,m, young, poisson));
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

void MPM::AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson)
{
	Vector3i maxx = Vector3i::Zero();
	for (int d=0; d<D; ++d)
	{
		maxx(d) = (l(d)/ratio)-1;
	}

	for (int k=0; k<=maxx(2); ++k)
    for (int j=0; j<=maxx(1); ++j)
    for (int i=0; i<=maxx(0); ++i)
    {
    	int type = -1;

    	if (D==1)
    	{
    		if (i==0 || i==maxx(0))	type = -2;
    	}
    	else if (D==2)
    	{
    		if (i==0 || i==maxx(0) || j==0 || j==maxx(1))	type = -2;		
    	}
    	else if (D==3)
    	{
    		if (i==0 || i==maxx(0) || j==0 || j==maxx(1) || k==0 || k==maxx(2))		type = -2;
    	}

    	Vector3d x = Vector3d::Zero();
    					x(0) = ratio*(i + 0.5)+x0(0);
    	if (D>1)		x(1) = ratio*(j + 0.5)+x0(1);
    	if (D>2)		x(2) = ratio*(k + 0.5)+x0(2);

    	AddParticle(type, x, m, young, poisson);
    }

    Lbp.resize(0);
    for (size_t p=0; p<Lp.size(); ++p)
    {
    	Lp[p]->Vol0 = 1.;
    	for (int d=0; d<D; ++d)
    	{
    		Lp[p]->Vol0 *= ratio;
    		if (Ntype==3)	Lp[p]->PSize0(d) = 0.5*ratio;
    		else 			Lp[p]->PSize0(d) = 0.;
    		Lp[p]->PSize(d) = Lp[p]->PSize0(d);
    	}
    	Lp[p]->Vol 	= Lp[p]->Vol0;

    	if (Lp[p]->Type==-2)
    	{
    		Lbp.push_back(Lp[p]);
    	}
    }
}

inline void MPM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};			//create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* tag_h5 	= new double[  Lp.size()];
	double* vel0_h5 = new double[  Lp.size()];
	double* pos_h5 	= new double[3*Lp.size()];
	double* vel_h5 	= new double[3*Lp.size()];

	for (size_t i=0; i<Lp.size(); ++i)
	{
        tag_h5[  i  ] 	= Lp[i]->Tag;
        // tag_h5[  i  ] 	= Lp[i]->ID;
        vel0_h5[ i  ] 	= Lp[i]->V.norm();
		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);
	}

	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_vel0	= new DataSet(file.createDataSet("Velocity_Magnitude", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
	dataset_vel0->write(vel0_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;
	delete dataset_tag;
	delete dataset_vel0;
	delete dataset_pos;
	delete dataset_vel;

	delete tag_h5;
	delete vel0_h5;
	delete pos_h5;
	delete vel_h5;

	file.close();

	string file_name_xmf = "particle_"+out.str()+".xmf";

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
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}