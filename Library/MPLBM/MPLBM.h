#include "../HEADER.h"
#include <MPM.h>
#include <LBM.h>

class MPLBM
{
public:
	MPLBM();
	~MPLBM();
	MPLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, int ntype, Vector3d dx);
	void Init(double rho0, Vector3d initV);
	void AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m);
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

inline MPLBM::MPLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, int ntype, Vector3d dx)
{
	DomLBM = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, nu);
	DomMPM = new MPM(ntype, nx, ny, nz, dx);

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

inline void MPLBM::AddBoxParticles(int tag, Vector3d& x0, Vector3d& l, double ratio, double m)
{
	DomMPM->AddBoxParticles(tag, x0, l, ratio, m);
}

inline void MPLBM::Solve(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			DomLBM->WriteFileH5(t,1);
			DomMPM->WriteFileH5(t);			
		}
		// cout << "1" << endl;
		DomLBM->CollideSRT();
		// cout << "2" << endl;
		DomLBM->ApplyWall();
		DomLBM->Stream();
		// cout << "3" << endl;
		// DomLBM->SetWall();
		DomLBM->CalRhoV();
		// cout << "4" << endl;
		DomMPM->ParticleToNode();
		// cout << "5" << endl;
		DomMPM->CalVOnNode();
		// cout << "6" << endl;
		// for (size_t n=0; n<DomMPM->LAn.size(); ++n)
		// {
		// 	size_t id = DomMPM->LAn[n];
	 //    	size_t i, j, k;
	 //    	DomMPM->FindIndex(id, i, j, k);
	 //    	// mpm node velocity 
	 //    	Vector3d vs = DomMPM->Ln[id]->V;
	 //    	// DomLBM->V[i][j][k] = vs;
	 //    	if (vs.norm()>0.01)
	 //    	{
	 //    		cout << "wrong vs" << endl;
	 //    		cout << vs.transpose() << endl;
	 //    		cout << DomMPM->Ln[id]->V.transpose() << endl;
	 //    		cout << i <<" " << j << " " << k << endl;
	 //    		abort();
	 //    	}
		// }
		for (size_t p=0; p<DomMPM->Lp.size(); ++p)
		{
			Vector3d vf = DomLBM->InterpolateV(DomMPM->Lp[p]->X);
			Vector3d vs = DomMPM->Lp[p]->V;
			Vector3d mv = DomLBM->Rho0*DomMPM->Lp[p]->Vol*(vs-vf);
			for (size_t l=0; l<DomMPM->Lp[p]->Lni.size(); ++l)
			{
				size_t id = DomMPM->Lp[p]->Lni[l];
		    	size_t i, j, k;
		    	DomMPM->FindIndex(id, i, j, k);
		    	double n = DomMPM->Lp[p]->LnN[l];
		    	DomLBM->ExForce[i][j][k] += n*mv;
			}
		}

		// for (size_t c=0; c<DomLBM->Lwall.size(); ++c)
		// {
		// 	size_t i = DomLBM->Lwall[c](0);
		// 	size_t j = DomLBM->Lwall[c](1);
		// 	size_t k = DomLBM->Lwall[c](2);
		// 	DomLBM->V[i][j][k].setZero();
		// }
		// cout << "7" << endl;
		DomMPM->NodeToParticle();

		// for (size_t n=0; n<DomMPM->LAn.size(); ++n)
		// {
		// 	size_t id = LAn[n];
	 //    	size_t i, j, k;
	 //    	DomMPM->FindIndex(id, i, j, k);

	 //    	DomMPM->Ln[id]->V

		// 	double vos = DomMPM->M[i][j][k]/Rhop;
		// 	double vof = 1.-vos;
		// 	if (vof>1. || vof <0.)
		// 	{
		// 		cout << "vof= " << vof << endl;
		// 		abort();
		// 	}
		// 	double re = vof*Dp*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]).norm()/DomLBM->Nu;
		// 	double vof2 = vof*vof;
		// 	double cd = 10.*vos/vof2 + vof2*(1.+1.5*sqrt(vos));
		// 	// cout << "cd= " << cd << endl;
		// 	if (re>0.01)	cd += 0.413*re/24./vof2*(1./vof+3.*vos*vof+8.4*pow(re, -0.343))/(1.+pow(10.,3.*vos)*pow(re,-0.5*(1.+4.*vos)));
		// 	// Vector3d drag = -18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]);
		// 	Vector3d drag = -3.*M_PI*(1.-vos)*DomLBM->Nu*DomLBM->Rho0*Dp*cd*(DomMPM->V[i][j][k]-DomLBM->V[i][j][k]);
		// 	// modify the soil momentum
		// 	if (vos>1.0e-12)
		// 	{
		// 		DomMPM->F [i][j][k] += drag;
		// 		DomMPM->Mv[i][j][k] += drag;
		// 	}
		// 	// if (DomMPM->Mv[i][j][k].norm()/DomMPM->M[i][j][k]>0.1)
		// 	// if (abs(drag(2))>0.15)
		// 	// {
		// 	// 	cout << "vos= " << vos << endl;
		// 	// 	cout << "vof= " << vof << endl;
		// 	// 	cout << "re= " << re << endl;
		// 	// 	cout << "cd= " << cd << endl;
		// 	// 	cout << "DomMPM->V[i][j][k]= " << DomMPM->V[i][j][k].transpose() << endl;
		// 	// 	cout << "DomLBM->V[i][j][k]= " << DomLBM->V[i][j][k].transpose() << endl;
		// 	// 	cout << "drag= " << drag.transpose() << endl;
		// 	// 	cout << "18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd= " << 18.*vos*vof*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd << endl;
		// 	// 	cout << "Dp= " << Dp << endl;
		// 	// 	cout << "DomLBM->Nu= " << DomLBM->Nu << endl;
		// 	// 	cout << "DomLBM->Rho0= " << DomLBM->Rho0 << endl;
		// 	// 	cout << "(DomMPM->V[i][j][k]-DomLBM->V[i][j][k])= " << (DomMPM->V[i][j][k]-DomLBM->V[i][j][k]).transpose() << endl;
		// 	// 	cout << "DomMPM->v= " << DomMPM->Mv[i][j][k].transpose()/DomMPM->M[i][j][k] << endl;
		// 	// 	cout << "DomMPM->M[i][j][k]= " << DomMPM->M[i][j][k] << endl;
		// 	// 	cout << "====================================================" << endl;
		// 	// 	abort();
		// 	// }

		// 	DomLBM->BodyForceLocal(i, j, k, -drag);
		// 	double vos = DomMPM->M[i][j][k]/Rhop;
		// 	double vof = 1.-vos;
		// 	// vof = 1.;
		// 	// volume average velocity
		// 	// Vector3d vm = (DomMPM->Mv[i][j][k] + vof*DomLBM->Rho[i][j][k]*DomLBM->V[i][j][k])/(DomMPM->M[i][j][k] + vof*DomLBM->Rho[i][j][k]);
		// 	Vector3d vm = vos*DomMPM->Mv[i][j][k]/DomMPM->M[i][j][k] + vof*DomLBM->V[i][j][k];
		// 	// modify the soil momentum
		// 	// Vector3d force = DomMPM->M[i][j][k]*vm - DomMPM->Mv[i][j][k];
		// 	// Vector3d force = -vof*DomLBM->Rho[i][j][k]*(vm-DomLBM->V[i][j][k]);
		// 	// DomMPM->F [i][j][k] += force;
		// 	// DomMPM->Mv[i][j][k] += force;
		// 	// modify the fluid momentum
		// 	DomLBM->V[i][j][k] = vm;
		// 	// DomLBM->BodyForceLocal(i, j, k, -force);
		// 	// cout << "DomMPM->M[i][j][k]= " << DomMPM->M[i][j][k] << endl;


		// 	// // volume average velocity
		// 	// Vector3d vm = (DomMPM->Mv[i][j][k] + vof*DomLBM->Rho[i][j][k]*DomLBM->V[i][j][k])/(DomMPM->M[i][j][k] + vof*DomLBM->Rho[i][j][k]);
		// 	// // modify the soil momentum
		// 	// DomMPM->F [i][j][k] += DomMPM->M[i][j][k]*vm - DomMPM->Mv[i][j][k];
		// 	// DomMPM->Mv[i][j][k] = DomMPM->M[i][j][k]*vm;
		// 	// // modify the fluid momentum
		// 	// DomLBM->V[i][j][k] = vm;
	 //    }
	}
}