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
		for (size_t p=0; p<DomMPM->Lp.size(); ++p)
		{
			Vector3d vf = DomLBM->InterpolateV(DomMPM->Lp[p]->X);	// Fluid velocity
			Vector3d vs = DomMPM->Lp[p]->V;							// Solid velocity
			Vector3d vsf = vs-sf;									// Relative velocity

			double n = DomMPM->M/Rhop/DomMPM->Vol;					// Porosity
			double n2 = n*n;
			double phi = 1.-n;
			double re = n*Dp*vsf.norm()/DomLBM->Nu;					// Reynolds number

			double cd = 10.*phi/n2 + n2*(1.+1.5*sqrt(phi));			// Drag coefficient
			cd += 0.03441666*re/n2*(1./n+3.*n*phi+8.4*pow(re, -0.343))/(1.+pow(10.,3.*phi)*pow(re,-0.5*(1.+4.*phi)));

			Vector3d drag = 18.*phi*n*DomLBM->Nu*DomLBM->Rho0/Dp/Dp*cd*vsf;	// volume averaged drag force
			drag *= DomMPM->Vol;

			Vector3d dragf = drag/n;

			for (size_t l=0; l<DomMPM->Lp[p]->Lni.size(); ++l)
			{
				size_t id = DomMPM->Lp[p]->Lni[l];
		    	size_t i, j, k;
		    	DomMPM->FindIndex(id, i, j, k);
		    	double ns = DomMPM->Lp[p]->LnN[l];
		    	DomLBM->ExForce[i][j][k] -= ns*dragf;
			}
		}
		DomMPM->NodeToParticle();
	}
}