#include "../HEADER.h"
#include <MPM.h>
#include <LBM.h>

class MPLBM
{
public:
	MPLBM();
	~MPLBM();
	MPLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, int ntype, int cmtype, Vector3d dx);
	void Init(double rho0, Vector3d initV);
	void AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson);
	void Solve(int tt, int ts);

	LBM* 							DomLBM;												// Domain of LBM
	MPM*							DomMPM;												// Domain of MPM

    int 							Nx;													// Domain size
    int 							Ny;
    int 							Nz;	

	int 							Nproc;												// Number of processors which used
    int 							D;

    double 							Rhop;												// Density of solid particle
    double 							Dp;													// Diameter of solid particle
};

inline MPLBM::MPLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, int ntype, int cmtype, Vector3d dx)
{
	DomLBM = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, nu);
	DomMPM = new MPM(ntype, cmtype, nx, ny, nz, dx);

	Nx = nx;
	Ny = ny;
	Nz = nz;

    D = DomLBM->D;
}

inline void MPLBM::Init(double rho0, Vector3d initV)
{
	DomLBM->Init(rho0, initV);
	DomMPM->Init();
}

inline void MPLBM::AddBoxParticles(Vector3d& x0, Vector3d& l, double ratio, double m, double young, double poisson)
{
	DomMPM->AddBoxParticles(x0, l, ratio, m, young, poisson);
}

inline void MPLBM::Solve(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			DomLBM->WriteFileH5(t);
			DomMPM->WriteFileH5(t);			
		}

		DomLBM->CollideMRT();
		DomLBM->Stream();
		// DomLBM->SetWall();
		DomLBM->CalRhoV();

		for (size_t c=0; c<DomMPM->LFn.size(); ++c)
		{
			int i = DomMPM->LFn[c][0];
			int j = DomMPM->LFn[c][1];
			int k = DomMPM->LFn[c][2];
			DomLBM->SBounceBack(i,j,k);
		}		

		DomMPM->ParticleToNode();
		if (t==0)	DomMPM->CalFOnNode(true);
		else 		DomMPM->CalFOnNode(false);

		for (size_t c=0; c<DomMPM->LFn.size(); ++c)
		{
			int i = DomMPM->LFn[c][0];
			int j = DomMPM->LFn[c][1];
			int k = DomMPM->LFn[c][2];
			int d = DomMPM->LFn[c][3];
			DomMPM->F[i][j][k](d) = 0.;
			DomMPM->Mv[i][j][k](d) = 0.;
		}

		for (size_t c=0; c<DomMPM->LAn.size(); ++c)
		{
			int i = DomMPM->LAn[c](0);
			int j = DomMPM->LAn[c](1);
			int k = DomMPM->LAn[c](2);

/*			double vos = DomMPM->M[i][j][k]/Rhop;
			double vof = 1.-vos;
			if (vof>1. || vof <0.)
			{
				cout << "vof= " << vof << endl;
				abort();
			}
			double re = vof*Dp*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]).norm()/DomLBM->Nu;
			double vof2 = vof*vof;
			double cd = 10.*vos/vof2 + vof2*(1.+1.5*sqrt(vos));
			// cout << "cd= " << cd << endl;
			if (re>0.01)	cd += 0.413*re/24./vof2*(1./vof+3.*vos*vof+8.4*pow(re, -0.343))/(1.+pow(10.,3.*vos)*pow(re,-0.5*(1.+4.*vos)));
			// Vector3d drag = -18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]);
			Vector3d drag = -3.*M_PI*(1.-vos)*DomLBM->Nu*DomLBM->Rho0*Dp*cd*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]);
			// modify the soil momentum
			if (vos>1.0e-12)
			{
				DomMPM->F [i][j][k] += drag;
				DomMPM->Mv[i][j][k] += drag;
			}
			// if (DomMPM->Mv[i][j][k].norm()/DomMPM->M[i][j][k]>0.1)
			// if (abs(drag(2))>0.15)
			// {
			// 	cout << "vos= " << vos << endl;
			// 	cout << "vof= " << vof << endl;
			// 	cout << "re= " << re << endl;
			// 	cout << "cd= " << cd << endl;
			// 	cout << "DomMPM->V[i][j][k]= " << DomMPM->V[i][j][k].transpose() << endl;
			// 	cout << "DomLBM->V[i][j][k]= " << DomLBM->V[i][j][k].transpose() << endl;
			// 	cout << "drag= " << drag.transpose() << endl;
			// 	cout << "18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd= " << 18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd << endl;
			// 	cout << "Dp= " << Dp << endl;
			// 	cout << "DomLBM->Nu= " << DomLBM->Nu << endl;
			// 	cout << "DomLBM->Rho0= " << DomLBM->Rho0 << endl;
			// 	cout << "(DomMPM->V[i][j][k]-DomLBM->V[i][j][k])= " << (DomMPM->V[i][j][k]-DomLBM->V[i][j][k]).transpose() << endl;
			// 	cout << "DomMPM->v= " << DomMPM->Mv[i][j][k].transpose()/DomMPM->M[i][j][k] << endl;
			// 	cout << "DomMPM->M[i][j][k]= " << DomMPM->M[i][j][k] << endl;
			// 	cout << "====================================================" << endl;
			// 	abort();
			// }

			DomLBM->BodyForceLocal(i, j, k, -drag);*/
			double vos = DomMPM->M[i][j][k]/Rhop;
			double vof = 1.-vos;
			// vof = 1.;
			// volume average velocity
			// Vector3d vm = (DomMPM->Mv[i][j][k] + vof*DomLBM->Rho[i][j][k]*DomLBM->V[i][j][k])/(DomMPM->M[i][j][k] + vof*DomLBM->Rho[i][j][k]);
			Vector3d vm = vos*DomMPM->Mv[i][j][k]/DomMPM->M[i][j][k] + vof*DomLBM->V[i][j][k];
			// modify the soil momentum
			// Vector3d force = DomMPM->M[i][j][k]*vm - DomMPM->Mv[i][j][k];
			// Vector3d force = -vof*DomLBM->Rho[i][j][k]*(vm-DomLBM->V[i][j][k]);
			// DomMPM->F [i][j][k] += force;
			// DomMPM->Mv[i][j][k] += force;
			// modify the fluid momentum
			DomLBM->V[i][j][k] = vm;
			// DomLBM->BodyForceLocal(i, j, k, -force);
			// cout << "DomMPM->M[i][j][k]= " << DomMPM->M[i][j][k] << endl;


			// // volume average velocity
			// Vector3d vm = (DomMPM->Mv[i][j][k] + vof*DomLBM->Rho[i][j][k]*DomLBM->V[i][j][k])/(DomMPM->M[i][j][k] + vof*DomLBM->Rho[i][j][k]);
			// // modify the soil momentum
			// DomMPM->F [i][j][k] += DomMPM->M[i][j][k]*vm - DomMPM->Mv[i][j][k];
			// DomMPM->Mv[i][j][k] = DomMPM->M[i][j][k]*vm;
			// // modify the fluid momentum
			// DomLBM->V[i][j][k] = vm;
	    }
		// DomLBM->Stream();
		// // // DomLBM->SetWall();
		// DomLBM->CalRhoV();
		DomMPM->CalVOnNode();
		DomMPM->NodeToParticle();
		if 		(DomMPM->CMType==0) 	DomMPM->CalStressOnParticleElastic();
		else if (DomMPM->CMType==1) 	DomMPM->CalStressOnParticleMohrCoulomb();
		else if (DomMPM->CMType==2) 	DomMPM->CalStressOnParticleNewtonian();
	}
}