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
#include <DEM_PARTICLE.h>

class DEM
{
public:
	DEM(int nx, int ny, int nz, int cmtype, double cr);
	~DEM();
	void Init();
	void AddSphere(int tag, double r, Vector3d& x, double rho);
	void AddCuboid(int tag, double lx, double ly, double lz, Vector3d& x, double rho);
	void AddDisk2D(int tag, double r, Vector3d& x, double rho);
	void AddNSpheres(int tag, int np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho);
	void Move();
	void ZeroForceTorque();
	void SetG(Vector3d& g);
	double EffectiveValue(double ai, double aj);											// Calculate effective values for contact force
	void RecordX();																			// Record position at Xb for check refilling LBM nodes
	void Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, bool& contacted);
	void UpdateFlag(DEM_PARTICLE* p0);
	void FindContact();
	void Contact();
	void LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void Solve(int tt, int ts, double dt, bool writefile);
	void LoadDEMFromH5( string fname);
	void WriteFileH5(int n);
	void WriteFileParticleInfo(int n);

	void (DEM::*ContactPara)(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);

	bool 							Periodic[3];

	int 							Nproc;
    int 							Nx;														// Mesh size for contact detection
    int 							Ny;
    int 							Nz;
    int 							CMType;

    int 							DomSize[3];

    double 							Dt;														// Time step
    double 							Cr;														// Coefficient of restitution
	double 							Beta;					
	double 							RatioGnt;												// Gt/Gn

    vector<size_t>***		 		Flag;													// Flag of lattice type

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < vector< size_t > >     Lc;                                                 	// List of potential contacted paricles' ID
	
	unordered_map<size_t, bool> 	CMap;													// Contact Map
	unordered_map<size_t, Vector3d> FMap;													// Friction Map
};

inline DEM::DEM(int nx, int ny, int nz, int cmtype, double cr)
{
	Nx = nx;
	Ny = ny;
	Nz = nz;

	DomSize[0] = Nx;
	DomSize[1] = Ny;
	DomSize[2] = Nz;

	Nproc = 1;

	Periodic[0] = true;
	Periodic[1] = true;
	Periodic[2] = true;

	Cr = cr;
	Beta = log(Cr)/sqrt(log(Cr)*log(Cr)+9.8696044010893586188344909998761);
	RatioGnt = 0.;

	CMType = cmtype;

	if (CMType==0)
	{
		ContactPara =& DEM::LinearContactPara;
		cout << "Using Linear contact model." << endl;

	}
	else if (CMType==1)
	{
		ContactPara =& DEM::HertzContactPara;
		cout << "Using Hertz Contact model." << endl;

	}
}

// Map a pair of integer to one integer key for hashing
// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
inline size_t Key(int i, int j)
{
	return (i+j+1)*(i+j)/2+j;
}

inline void DEM::Init()
{
    cout << "================ Start init. ================" << endl;

    Dt = 1.;

	Flag = new vector<size_t>** [Nx+1];

	for (int i=0; i<=Nx; ++i)
	{
		Flag[i] = new vector<size_t>* [Ny+1];

		for (int j=0; j<=Ny; ++j)
		{
			Flag[i][j]	= new vector<size_t> [Nz+1];

			for (int k=0; k<=Nz; ++k)
			{
				Flag[i][j][k].resize(0);
			}
		}
	}
	// Mark flag and add particles for the surounding box to handle walls and periodic BC
	// For x axis
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		Flag[0 ][j][k].push_back(0);
		Flag[Nx][j][k].push_back(1);
	}
	Vector3d x0 (0., 0., 0.);
	Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
	Lp[0]->ID = 0;
	Lp[0]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
	Lp[1]->ID = 1;
	Lp[1]->Type = -1;
	// For y axis
	for (int i=0; i<=Nx; ++i)
	for (int k=0; k<=Nz; ++k)
	{
		Flag[i][0 ][k].push_back(2);
		Flag[i][Ny][k].push_back(3);
	}		
	Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
	Lp[2]->ID = 2;
	Lp[2]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
	Lp[3]->ID = 3;
	Lp[3]->Type = -1;
	// For z axis
	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	{
		Flag[i][j][0 ].push_back(4);
		Flag[i][j][Nz].push_back(5);
	}
	Lp.push_back(new DEM_PARTICLE(-5, x0, 0.));
	Lp[4]->ID = 4;
	Lp[4]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-6, x0, 0.));
	Lp[5]->ID = 5;
	Lp[5]->Type = -1;
	cout << "================ Finish init. ================" << endl;
}

inline void DEM::AddSphere(int tag, double r, Vector3d& x, double rho)
{
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
	Lp[Lp.size()-1]->SetSphere(r);
	UpdateFlag(Lp[Lp.size()-1]);
}

inline void DEM::AddCuboid(int tag, double lx, double ly, double lz, Vector3d& x, double rho)
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		if (tag==Lp[i]->Tag)
		{
			cout << "\033[1;31mError: Cannot add sphere. The tag is existed.\033[0m\n";		
			exit(0);
		}
	}
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
	cout << "SetCuboid" << endl;
	Lp[Lp.size()-1]->SetCuboid(lx, ly, lz);
	cout << "SetCuboid Finish" << endl;
}

inline void DEM::AddDisk2D(int tag, double r, Vector3d& x, double rho)
{
	if (x(2)!=0.)
	{
		cout << "\033[1;31mError: Disk only works on 2D, check position again.\033[0m\n";		
		exit(0);		
	}

	for (size_t i=0; i<Lp.size(); ++i)
	{
		if (tag==Lp[i]->Tag)
		{
			cout << "\033[1;31mError: Cannot add sphere. The tag is existed.\033[0m\n";		
			exit(0);
		}
	}
	Lp.push_back(new DISK2D(tag, r, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
}

inline void DEM::AddNSpheres(int tag, int np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho)
{
    srand (time(NULL));
    int count = 0;

    Vector3d l = x1-x0;
    for (int n=0; n<1.0e200; n++)
    {
        Vector3d pos(0.,0.,0.);
        pos(0) = (l(0)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        pos(1) = (l(1)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        pos(2) = (l(2)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);

        bool boolContact = false;

        if (Lp.size()!=0)
        {
            // #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t p=0; p<Lp.size(); ++p)
            {
                Vector3d xpq = Lp[p]->X - pos;
                double pqDis = xpq.norm() - r - Lp[p]->R;
                if (pqDis<surDis)    boolContact = true;           
            }
        }
        if (boolContact==false)   
        {
            pos += x0;
            AddSphere(tag, r, pos, rho);
            count++;
            cout << "Already add " << count << " particles" << endl;
        }

        if (count>=np)   break;
    }
}

inline void DEM::Move()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		DEM_PARTICLE* p0 = Lp[i];
		if 		(p0->X(0)>Nx)	p0->X(0) = p0->X(0)-Nx-1;
		else if (p0->X(0)<0.)	p0->X(0) = p0->X(0)+Nx+1;
		if 		(p0->X(1)>Ny)	p0->X(1) = p0->X(1)-Ny-1;
		else if (p0->X(1)<0.)	p0->X(1) = p0->X(1)+Ny+1;
		if 		(p0->X(2)>Nz)	p0->X(2) = p0->X(2)-Nz-1;
		else if (p0->X(2)<0.)	p0->X(2) = p0->X(2)+Nz+1;
		Lp[i]->VelocityVerlet(Dt);
	}
}

inline void DEM::ZeroForceTorque()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		Lp[i]->ZeroForceTorque();
	}
}

inline void DEM::SetG(Vector3d& g)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		Lp[i]->SetG(g);
	}
}

inline void DEM::RecordX()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		Lp[i]->Xb = Lp[i]->X;
	}
}

void DEM::LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt)
{
	double me 	= 2.*EffectiveValue(pi->M, pj->M);

	// For collision with wall
	if (pi->Type==-1)
	{
		me = pj->M;
	}

	kn 	= 2.*EffectiveValue(pi->Kn, pj->Kn);
	gn 	= -2.*Beta*sqrt(kn*me);
	gn 	= 2.*EffectiveValue(pi->Gn, pj->Gn);
	kt 	= 2.*EffectiveValue(pi->Kt, pj->Kt);
	gt 	= RatioGnt*gn;	
}

void DEM::HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt)
{
	double re 	= EffectiveValue(pi->R, pj->R);
	double me 	= EffectiveValue(pi->M, pj->M);
	double ee 	= 1./((1.-pi->Poisson*pi->Poisson)/pi->Young + (1.-pj->Poisson*pj->Poisson)/pj->Young);
	double ge 	= 1./(2.*(2.-pi->Poisson)*(1.+pi->Poisson)/pi->Young + 2.*(2.-pj->Poisson)*(1.+pj->Poisson)/pj->Young);

	// For collision with wall
	if (pi->Type==-1)
	{
		me = pj->M;
		re = pj->R;
	}

	double sn 	= 2.*ee*sqrt(re*delta);
	double st 	= 8.*ge*sqrt(re*delta);
	kn 	= 1.3333333333333*ee*sqrt(re*delta);
	gn 	= -1.825741858351*Beta*sqrt(sn*me);
	kt 	= 8.*ge*sqrt(re*delta);
	gt 	= -1.825741858351*Beta*sqrt(st*me);	
	gn = -log(Cr)*me;
	gt = 0.;
}

// Contact force model for spheres
inline void DEM::Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, bool& contacted)
{
	contacted = false;
	Vector3d Xi = pi->X;
	Vector3d Xj = pj->X;

	if (pi->Type==-1)
	{
		Xi = Xj;
		int axis = pi->ID/2;
		int dirc = pi->ID-2*axis;
		Xi(axis) = dirc*DomSize[axis];
	}

	// Normal direction (pj pinnts to pi)
	Vector3d n = Xi-Xj;
	// Overlapping distance
	double delta = pi->R+pj->R-n.norm();

	if (pi->crossingFlag || pj->crossingFlag)
	{
		vector<Vector3d> LXimir;
		LXimir.push_back(Xi);
		vector<Vector3d> LXjmir;
		LXimir.push_back(Xj);

		Vector3d Ximir = Xi;
		Vector3d Xjmir = Xj;

		for (int d=0; d<3; ++d)
		{
			if (pi->crossing[d])
			{
				if (Xi(d)+pi->R>DomSize[d])		Ximir(d) = Xi(d)-DomSize[d]-1;
				else if (Xi(d)-pi->R<0.)		Ximir(d) = Xi(d)+DomSize[d]+1;

			}
			if (pj->crossing[d])
			{
				if (Xj(d)+pj->R>DomSize[d])		Xjmir(d) = Xj(d)-DomSize[d]-1;
				else if (Xj(d)-pj->R<0.)		Xjmir(d) = Xj(d)+DomSize[d]+1;
			}
		}

		double Xit[2][3] = {{Xi(0), Xi(1), Xi(2)}, {Ximir(0), Ximir(1), Ximir(2)}};
		double Xjt[2][3] = {{Xj(0), Xj(1), Xj(2)}, {Xjmir(0), Xjmir(1), Xjmir(2)}};

		int loopOrder[7][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,0,1}, {0,1,1}, {1,1,1}};

		vector<Vector3d> Lxi, Lxj;
		Lxi.push_back(Xi);
		Lxj.push_back(Xj);

		for (int m=0; m<7; ++m)
		{
			int i = loopOrder[m][0];
			int j = loopOrder[m][1];
			int k = loopOrder[m][2];
			Vector3d Xim (Xit[i][0], Xit[j][1], Xit[k][2]);
			if (Xim!=Xi)	Lxi.push_back(Xim);
			Vector3d Xjm (Xjt[i][0], Xjt[j][1], Xjt[k][2]);
			if (Xjm!=Xj)	Lxj.push_back(Xjm);
		}

		for (size_t i=0; i<Lxi.size(); ++i)
		for (size_t j=0; j<Lxj.size(); ++j)
		{
			double deltat = pi->R+pj->R-(Lxi[i]-Lxj[j]).norm();
			if (deltat>0.)
			{
				delta = deltat;
				n = Lxi[i]-Lxj[j];
				break;
			}
		}
	}

	if (delta>0.)
	{
		if (delta>0.01*(pi->R+pj->R))
		{
			cout << "max delta= " << 0.01*(pi->R+pj->R) << endl;
			cout << "xi= " << Xi << endl;
			cout << "xj= " << Xj << endl;
			cout << "delta= " << delta << endl;
			cout << pi->ID << endl;
			cout << pj->ID << endl;
			cout << Xi.transpose() << endl;
			cout << Xj.transpose() << endl;
			cout << pi->Xmir.transpose() << endl;
			cout << pj->Xmir.transpose() << endl;
			abort();			
		}
		// report contact for updating friction map 
		contacted = true;
		n.normalize();

		double kn, gn, kt, gt;
		(this->*ContactPara)(pi, pj, delta, kn, gn, kt, gt);
 		// Relative velocity in normal direction
		Vector3d vn = (pj->V-pi->V).dot(n)*n;
		// Normal contact force
		// Vector3d fn= kn*delta*n + 2.*gn*sqrt(0.5*me*kn)*vn;
		Vector3d fn= kn*delta*n + gn*vn;
		// Relative velocity at the contact point
		Vector3d vij = pi->V-pj->V+(pi->R-0.5*delta)*n.cross(pi->W)+(pj->R-0.5*delta)*n.cross(pj->W);	// eq.10
		// Relative tangential velocity at the contact point
		Vector3d vt = vij-n.dot(vij)*n;						// eq.9
		// Update tangential spring
		Vector3d xi0 = vt*Dt;
		size_t key = Key(pi->ID, pj->ID);
		if (FMap.count(key)>0)	xi0 += FMap[key];			// eq.19
		// Project to current tangential plane
		xi = xi0 - n.dot(xi0)*n;							// eq.17
		// Pick the larger friction coefficient
		double mu_s = 0.;
		if (pi->Mu_s>1.0e-12 && pj->Mu_s>1.0e-12)	mu_s = max(pi->Mu_s, pj->Mu_s);
		// double mu_d = 0.;
		// if (pi->Mu_d>1.0e-12 && pj->Mu_d>1.0e-12)	mu_d = max(pi->Mu_d, pj->Mu_d);	
		// Static tangential force
		double fts = mu_s*fn.norm();
		// Tangential force
		Vector3d ft = -kt*xi-gt*vt;							// eq.18
		if (ft.norm()>fts)
		{
			Vector3d t = ft.normalized();
			ft = fts*t;										// eq.21
			xi = -(fts*t+gt*vt)/kt;							// eq.20
		}
		// Torque with normalized arm
		Vector3d arm = -n.cross(ft);
		// Only works for sphere
		#pragma omp critical
		{
			pi->Fc += fn+ft;
			pj->Fc -= fn+ft;

			pi->Tc += (pi->R-0.5*delta)*arm;
			pj->Tc += (pj->R-0.5*delta)*arm;
		}
	}
}

inline double DEM::EffectiveValue(double ai, double aj)
{
	double a = 0.;
	if (ai>1.0e-12 && aj>1.0e-12)	a = ai*aj/(ai+aj);
	return (a);
}

inline void DEM::UpdateFlag(DEM_PARTICLE* p0)
{
	for (int i=p0->MinX; i<=p0->MaxX; ++i)
	for (int j=p0->MinY; j<=p0->MaxY; ++j)
	for (int k=p0->MinZ; k<=p0->MaxZ; ++k)
	{
		Vector3d ind (i,j,k);
		// Distance from the cell centre to particle surface
		double dis = 0.;
		if (p0->Type==1)	dis = (ind-p0->X).norm()-p0->R-0.87;
		if (dis<0)
		{
			int ic = (i+Nx+1)%(Nx+1);
			int jc = (j+Ny+1)%(Ny+1);
			int kc = (k+Nz+1)%(Nz+1);

			size_t p = p0->ID;
			if (Flag[ic][jc][kc].size()>0)
			{
				for (size_t m=0; m<Flag[ic][jc][kc].size(); ++m)
				{
					size_t q = Flag[ic][jc][kc][m];
					if (q<6 && Periodic[q/2])
					{
						p0->crossing[q/2] = true;
						p0->crossingFlag = true;
					}
					else
					{
						size_t min0 = min(p,q);
						size_t max0 = max(p,q);
						size_t key = Key(min0,max0);
						if (!CMap[key])
						{
							Lc.push_back({min0, max0});
							CMap[key] = true;
						}
					}		
				}
			}
			Flag[ic][jc][kc].push_back(p);
		}
	}	
}

inline void DEM::FindContact()
{
	CMap.clear();
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];

		for (int i=p0->MinX; i<=p0->MaxX; ++i)
		for (int j=p0->MinY; j<=p0->MaxY; ++j)
		for (int k=p0->MinZ; k<=p0->MaxZ; ++k)
		{
			int ic = (i+Nx+1)%(Nx+1);
			int jc = (j+Ny+1)%(Ny+1);
			int kc = (k+Nz+1)%(Nz+1);

			if (Flag[ic][jc][kc].size()>0)
			{
				vector<size_t> flagtmp;
				flagtmp.resize(0);
				for (size_t m=0; m<Flag[ic][jc][kc].size(); ++m)
				{
					size_t ind = Flag[ic][jc][kc][m];
					if (Lp[ind]->fixed || ind<6)
					{
						flagtmp.push_back(ind);
					}
				}
				Flag[ic][jc][kc] = flagtmp;
			}
		}
	}
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];
		if (!p0->fixed)	UpdateFlag(p0);
	}
}

inline void DEM::Contact()
{
    if (Lc.size()>0.)
    {
    	unordered_map<size_t, Vector3d> fmap;
    	// #pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t l=0; l<Lc.size(); ++l)
		{
			int i = Lc[l][0];
			int j = Lc[l][1];

			bool contacted = false;
			Vector3d xi (0.,0.,0.);
			Contact2P(Lp[i], Lp[j], xi, contacted);
			if (contacted)	fmap[Key(i,j)] = xi;
		}
		FMap = fmap;
    }
	Lc.clear();
}

inline void DEM::Solve(int tt, int ts, double dt, bool writefile)
{
	Dt = dt;
	bool show = 1;
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step ============ " << t << endl;
			if (writefile)	WriteFileH5(t);
		}

		if (show)	cout << "Time Step ============ " << t << endl;
		clock_t t_start = std::clock();
		FindContact();
		clock_t t_end = std::clock();
		if (show)	cout << "FindContact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::clock();
		Contact();
		t_end = std::clock();
		if (show)	cout << "Contact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (t==0)
		{
			for (size_t p=0; p<Lp.size(); ++p)
			{
				Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
				Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
			}
		}
		t_start = std::clock();
		Move();
		t_end = std::clock();
		if (show)	cout << "Move time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::clock();
		// WriteFileParticleInfo(t);
		ZeroForceTorque();
		t_end = std::clock();
		if (show)	cout << "ZeroForceTorque time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	}
}

void DEM::LoadDEMFromH5( string fname)
{
	cout << "========= Start loading DEM particles from " << fname << "==============" << endl;
	H5std_string FILE_NAME( fname );
	H5std_string DATASET_NAME_POS( "Position" );
	H5File file_pos( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_pos = file_pos.openDataSet( DATASET_NAME_POS );
	DataSpace dataspace_pos = dataset_pos.getSpace();
    hsize_t dims_pos[2];
    dataspace_pos.getSimpleExtentDims( dims_pos, NULL);
    hsize_t dimsm_pos = dims_pos[0];

	H5std_string DATASET_NAME_R( "Radius" );
	H5File file_r( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_r = file_r.openDataSet( DATASET_NAME_R );
	DataSpace dataspace_r = dataset_r.getSpace();
    hsize_t dims_r[2];
    dataspace_r.getSimpleExtentDims( dims_r, NULL);
    hsize_t dimsm_r = dims_r[0];

	H5std_string DATASET_NAME_RHO( "Rho" );
	H5File file_rho( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_rho = file_rho.openDataSet( DATASET_NAME_RHO );
	DataSpace dataspace_rho = dataset_r.getSpace();
    hsize_t dims_rho[2];
    dataspace_rho.getSimpleExtentDims( dims_rho, NULL);
    hsize_t dimsm_rho = dims_rho[0];

	H5std_string DATASET_NAME_TAG( "Tag" );
	H5File file_tag( FILE_NAME, H5F_ACC_RDONLY );
	DataSet dataset_tag = file_tag.openDataSet( DATASET_NAME_TAG );
	DataSpace dataspace_tag = dataset_tag.getSpace();
    hsize_t dims_tag[2];
    dataspace_tag.getSimpleExtentDims( dims_tag, NULL);
    hsize_t dimsm_tag = dims_tag[0];

    double data_pos[dimsm_pos];
    dataset_pos.read( data_pos, PredType::NATIVE_DOUBLE, dataspace_pos, dataspace_pos );

    double data_r[dimsm_r];
    dataset_r.read( data_r, PredType::NATIVE_DOUBLE, dataspace_r, dataspace_r );

    double data_rho[dimsm_rho];
    dataset_rho.read( data_rho, PredType::NATIVE_DOUBLE, dataspace_rho, dataspace_rho );

    double data_tag[dimsm_tag];
    dataset_tag.read( data_tag, PredType::NATIVE_DOUBLE, dataspace_tag, dataspace_tag );

    int np = dimsm_pos/3;
    for (int i=6; i<np; ++i)
    {
    	Vector3d pos (data_pos[3*i], data_pos[3*i+1], data_pos[3*i+2]);
    	double r = data_r[i];
    	int tag = (int) data_tag[i];
    	double rho = data_rho[i];
    	AddSphere(tag, r, pos, rho);
    }
    cout << "========= Loaded "<< Lp.size()-6<< " DEM particles from " << fname << "==============" << endl;
}

inline void DEM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "DEM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};			//create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.

	hsize_t n_points = 8*Lp.size();
	hsize_t n_faces  = 6*Lp.size();
	hsize_t	dims_points [1] = {3*n_points};			//create data space.
	hsize_t	dims_faces  [1] = {5*n_faces};			//create data space.
	hsize_t	dims_fscalar[1] = {n_faces};			//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	int rank_points = sizeof(dims_points ) / sizeof(hsize_t);
	int rank_faces  = sizeof(dims_faces  ) / sizeof(hsize_t);
	int rank_fscalar= sizeof(dims_fscalar) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	DataSpace	*space_points = new DataSpace(rank_points, dims_points);
	DataSpace	*space_faces  = new DataSpace(rank_faces , dims_faces );
	DataSpace	*space_fscalar= new DataSpace(rank_fscalar,dims_fscalar);

	double* r_h5 	= new double[  Lp.size()];
	double* rho_h5 	= new double[  Lp.size()];
	double* tag_h5 	= new double[  Lp.size()];
	double* pos_h5 	= new double[3*Lp.size()];
	double* vel_h5 	= new double[3*Lp.size()];
	double* agv_h5	= new double[3*Lp.size()];

	double* poi_h5	= new double[3*n_points];
	int* 	fac_h5	= new int   [5*n_faces];

	double* fv_h5 	= new double[n_faces];

	int count_p = 0;
	int count_f = 0;
	int count_fv = 0;

	for (size_t i=0; i<Lp.size(); ++i)
	{
		Vector3d agv = Lp[i]->Q._transformVector(Lp[i]->W);
        r_h5  [  i  ] 	= Lp[i]->R;
        rho_h5[  i  ] 	= Lp[i]->Rho;
        tag_h5[  i  ] 	= Lp[i]->Tag;
		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);
		agv_h5[3*i  ] 	= agv(0);
		agv_h5[3*i+1] 	= agv(1);
		agv_h5[3*i+2] 	= agv(2);

		for (size_t j=0; j<Lp[i]->P.size(); ++j)
		{
			poi_h5[count_p  ] = Lp[i]->P[j](0);
			poi_h5[count_p+1] = Lp[i]->P[j](1);
			poi_h5[count_p+2] = Lp[i]->P[j](2);
			count_p += 3;
		}
		// cout << "Lp[i]->Faces.size()= " << Lp[i]->Faces.size() << endl;
		for (size_t k=0; k<Lp[i]->Faces.size(); ++k)
		{
			fac_h5[count_f  ] = 5;
			fac_h5[count_f+1] = Lp[i]->Faces[k](0);
			fac_h5[count_f+2] = Lp[i]->Faces[k](1);
			fac_h5[count_f+3] = Lp[i]->Faces[k](2);
			fac_h5[count_f+4] = Lp[i]->Faces[k](3);
			count_f += 5;
			fv_h5[count_fv] = Lp[i]->V.norm();
			count_fv++;
		}
	}

	DataSet	*dataset_r 		= new DataSet(file.createDataSet("Radius", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_rho	= new DataSet(file.createDataSet("Rho", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_agv	= new DataSet(file.createDataSet("AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector));

    DataSet	*dataset_poi	= new DataSet(file.createDataSet("Points", PredType::NATIVE_DOUBLE, *space_points));
    DataSet	*dataset_fac	= new DataSet(file.createDataSet("Faces", PredType::NATIVE_INT, *space_faces));
    DataSet	*dataset_fv		= new DataSet(file.createDataSet("FaceVelocity", PredType::NATIVE_DOUBLE, *space_fscalar));

	dataset_r->write(r_h5, PredType::NATIVE_DOUBLE);
	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
	dataset_agv->write(agv_h5, PredType::NATIVE_DOUBLE);

	dataset_poi->write(poi_h5, PredType::NATIVE_DOUBLE);
	dataset_fac->write(fac_h5, PredType::NATIVE_INT);
	dataset_fv->write(fv_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;

	delete dataset_r;
	delete dataset_rho;
	delete dataset_tag;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_agv;

	delete dataset_poi;
	delete dataset_fac;
	delete dataset_fv;

	delete r_h5;
	delete rho_h5;
	delete tag_h5;
	delete pos_h5;
	delete vel_h5;
	delete agv_h5;

	delete poi_h5;
	delete fac_h5;
	delete fv_h5;

	file.close();

	string file_name_xmf = "DEM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";

    if (n_faces>0)
    {
	    oss << "   <Grid Name=\"DEM_FACES\">\n";
	    oss << "     <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << n_faces << "\">\n";
	    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << 5*n_faces << "\">\n";
	    oss << "        " << file_name_h5 <<":/Faces \n";
	    oss << "       </DataItem>\n";
	    oss << "     </Topology>\n";
	    oss << "     <Geometry GeometryType=\"XYZ\">\n";
	    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_points << " 3\" >\n";
	    oss << "        " << file_name_h5 <<":/Points \n";
	    oss << "       </DataItem>\n";
	    oss << "     </Geometry>\n";
	    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
	    oss << "       <DataItem Dimensions=\"" << n_faces << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
	    oss << "        " << file_name_h5 <<":/FaceVelocity\n";
	    oss << "       </DataItem>\n";
	    oss << "     </Attribute>\n";
	    oss << "   </Grid>\n";
    }

    oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Lp.size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Lp.size() << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Rho \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
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
    oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/AngularVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}

inline void DEM::WriteFileParticleInfo(int n)
{
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<Lp.size(); ++p)
    {
    	if (p==6 || p==6)
    	{
	        ostringstream info;
	        info << "particleInfo_" << p << ".res";
	        ofstream log(info.str().c_str(), ios_base::out | ios_base::app);
	        if (n==0)   log << "\"Time\"     \"X\"     \"Y\"     \"Z\"     \"U\"     \"V\"     \"W\"      \"Wx\"     \"Wy\"     \"Wz\"      \"Fhx\"     \"Fhy\"     \"Fhz\"      \"Thx\"     \"Thy\"     \"Thz\"     \"Fcx\"     \"Fcy\"     \"Fcz\"      \"Tcx\"     \"Tcy\"     \"Tcz\"\n";
	        else;
	        log << setprecision(9) << fixed << (double)n << "     " << Lp[p]->X(0) << "     " << Lp[p]->X(1) << "     " << Lp[p]->X(2) << "     " << Lp[p]->V(0) << "     " << Lp[p]->V(1) << "     " << Lp[p]->V(2) << "     " << Lp[p]->W(0) << "     " << Lp[p]->W(1) << "     " << Lp[p]->W(2) << "     " << Lp[p]->Fh(0) << "     " << Lp[p]->Fh(1) << "     " << Lp[p]->Fh(2)<< "     " << Lp[p]->Th(0) << "     " << Lp[p]->Th(1) << "     " << Lp[p]->Th(2) << "     " << Lp[p]->Fc(0) << "     " << Lp[p]->Fc(1) << "     " << Lp[p]->Fc(2)<< "     " << Lp[p]->Tc(0) << "     " << Lp[p]->Tc(1)<< "     " << Lp[p]->Tc(2) << "\n";

    	}
    }
}

// inline void DEM::WriteFileParticleInfo(int n)
// {
//     // #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (size_t p=6; p<Lp.size(); ++p)
//     {
//         ostringstream info;
//         info << "particleInfo_" << p << ".res";
//         ofstream log(info.str().c_str(), ios_base::out | ios_base::app);
//         if (n==0)   log << "\"Time\"     \"X\"     \"Y\"     \"Z\"     \"U\"     \"V\"     \"W\"      \"Wx\"     \"Wy\"     \"Wz\"      \"Fhx\"     \"Fhy\"     \"Fhz\"      \"Thx\"     \"Thy\"     \"Thz\"     \"Fcx\"     \"Fcy\"     \"Fcz\"      \"Tcx\"     \"Tcy\"     \"Tcz\"\n";
//         else;
//         log << setprecision(9) << fixed << (double)n << "     " << Lp[p]->X(0) << "     " << Lp[p]->X(1) << "     " << Lp[p]->X(2) << "     " << Lp[p]->V(0) << "     " << Lp[p]->V(1) << "     " << Lp[p]->V(2) << "     " << Lp[p]->W(0) << "     " << Lp[p]->W(1) << "     " << Lp[p]->W(2) << "     " << Lp[p]->Fh(0) << "     " << Lp[p]->Fh(1) << "     " << Lp[p]->Fh(2)<< "     " << Lp[p]->Th(0) << "     " << Lp[p]->Th(1) << "     " << Lp[p]->Th(2) << "     " << Lp[p]->Fc(0) << "     " << Lp[p]->Fc(1) << "     " << Lp[p]->Fc(2)<< "     " << Lp[p]->Tc(0) << "     " << Lp[p]->Tc(1)<< "     " << Lp[p]->Tc(2) << "\n";
//     }
// }