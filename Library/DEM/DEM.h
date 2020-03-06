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
#include <GJK.h>
// #include <2D_PDEM_FUNCTIONS.h>

class DEM
{
public:
	DEM(int nx, int ny, int nz, string cmtype, string dmtype, double cr);
	~DEM();
	void Init();
	void AddSphere(int tag, double r, Vector3d& x, double rho);
	void AddCuboid(int tag, double lx, double ly, double lz, Vector3d& x, double rho);
	void AddTetrahedron(int tag, vector<Vector3d> ver, double rho);
	void AddDisk2D(int tag, double r, Vector3d& x, double rho);
	void AddNSpheres(int tag, int np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho);
	void Add2DPolynomialParticle(int tag, VectorXd& coef, Vector3d& x, double rho);
	void Move();
	void ZeroForceTorque(bool h, bool c);
	void SetG(Vector3d& g);
	double EffectiveValue(double ai, double aj);											// Calculate effective values for contact force
	void RecordX();																			// Record position at Xb for check refilling LBM nodes
	void Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, Vector3d& xir, bool& contacted);
	void Friction(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kt, double gt, Vector3d& n, Vector3d& fn, Vector3d& xi, Vector3d& ft);
	void RollingResistance(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kr, double gr, Vector3d& n, Vector3d& fn, Vector3d& xir, Vector3d& armr);
	void UpdateFlag(DEM_PARTICLE* p0);
	void UpdateXmir(DEM_PARTICLE* p0);
	void UpdateXmirGlobal();
	void FindContact();
	void FindContactBasedOnNode(bool first);
	void Contact(bool writeFc, int n);
	void LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void LinearDampingPara0(double& kn, double& me, double& gn, double& gt);
	void HertzDampingPara0(double& kn, double& me, double& gn, double& gt);
	void DampingParaDoNothing(double& kn, double& me, double& gn, double& gt);
	void SetLubrication(double hn, double viscosity);
	void Solve(int tt, int ts, double dt, bool writefile);
	void DeleteParticles();
	void LoadDEMFromH5( string fname, double scale, double rhos);
	void WriteFileH5(int n);
	void WriteContactForceFileH5(int n);
	void WriteFileParticleInfo(int n);

	void (DEM::*ContactPara)(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void (DEM::*DampingPara)(double& kn, double& me, double& gn, double& gt);

	bool 							Periodic[3];

	size_t 							Nproc;
    size_t 							D;														// Dimension
    int 							Nx;														// Mesh size for contact detection
    int 							Ny;
    int 							Nz;
    string 							CMType;
    string 							DMType;

    int 							DomSize[3];
    size_t 							Np;														// Total number of points in the domain
    size_t 							Nf;														// Total number of faces in the domain
    double 							Dt;														// Time step
    double 							Cr;														// Coefficient of restitution
	double 							Beta;					
	double 							RatioGnt;												// Gt/Gn
	double 							RatioKnr;												// Kr/Kn
	double 							RatioGnr;												// Gr/Gn
	double 							Hn;														// Lubrication cut off distance for fluid
	double 							Viscosity;												// Dynamic viscosity fluid

	vector<Vector3i> 				Ne;														// Relative location of neighbor cells

    vector<size_t>***		 		Flag;													// Flag of lattice type
    vector<size_t>***		 		Flagt;													// Flag of lattice type

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < DEM_PARTICLE* >		Lg;														// List of groups
	vector < vector< size_t > >     Lc;                                                 	// List of potential contacted paricles' ID
	
	unordered_map<size_t, bool> 	CMap;													// Contact Map
	unordered_map<size_t, Vector3d> FMap;													// Friction Map
	unordered_map<size_t, Vector3d> RMap;													// Rolling resistance Map

	double 							FsTable[10][10];										// Friction coefficient (static) table
	double 							FdTable[10][10];										// Friction coefficient (dynamic) table
	double 							RsTable[10][10];										// Rolling friction coefficient (static) table
	double 							RdTable[10][10];										// Rolling Friction coefficient (dynamic) table
};

inline DEM::DEM(int nx, int ny, int nz, string cmtype, string dmtype, double cr)
{
	Nx = nx;
	Ny = ny;
	Nz = nz;
	Np = 0;
	Nf = 0;

	DomSize[0] = Nx;
	DomSize[1] = Ny;
	DomSize[2] = Nz;
	D = 3;
	if (Nz==0)
	{
		D = 2;
		if (Ny==0)	D = 1;
	}

	if (D==3)
	{
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };		
	}
	else if (D==2)
	{
		Ne  = {     { 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
	}

	Nproc = 1;

	Periodic[0] = true;
	Periodic[1] = true;
	Periodic[2] = true;

	Cr = cr;
	Beta = log(Cr)/sqrt(log(Cr)*log(Cr)+9.8696044010893586188344909998761);
	RatioGnt = 0.;

	Hn = 0.;
	Viscosity = 0.;

	CMType = cmtype;
	DMType = dmtype;

	if (CMType=="LINEAR")
	{
		ContactPara =& DEM::LinearContactPara;
		cout << "Using Linear contact model." << endl;
		if (DMType=="DEFULT")
		{
			DampingPara =& DEM::DampingParaDoNothing;
			cout << "Using defult Linear damping model." << endl;
		}
		else if (DMType=="CR")
		{
			DampingPara =& DEM::HertzDampingPara0;
			cout << "Using Coefficient of restitution for Linear damping model." << endl;
		}
		else
		{
			cout << "Undefined damping model!" << endl;
			abort();
		}
	}
	else if (CMType=="HERTZ")
	{
		ContactPara =& DEM::HertzContactPara;
		cout << "Using Hertz Contact model." << endl;
		if (DMType=="DEFULT")
		{
			DampingPara =& DEM::DampingParaDoNothing;
			cout << "Using defult Hertz damping model." << endl;
		}
		else if (DMType=="LOG")
		{
			DampingPara =& DEM::HertzDampingPara0;
			cout << "Using eq: log(Cr)*me for Hertz damping model." << endl;
		}
		else
		{
			cout << "Undefined damping model!" << endl;
			abort();
		}
	}
	Lp.resize(0);
	Lg.resize(0);
	// Set friction coefficient to zero
	for (size_t i=0; i<10; ++i)
	for (size_t j=0; j<10; ++j)
	{
		FsTable[i][j] = 0.;
		FdTable[i][j] = 0.;
		RsTable[i][j] = 0.;
		RdTable[i][j] = 0.;
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
	Flagt = new vector<size_t>** [Nx+1];

	for (int i=0; i<=Nx; ++i)
	{
		Flag[i] = new vector<size_t>* [Ny+1];
		Flagt[i] = new vector<size_t>* [Ny+1];

		for (int j=0; j<=Ny; ++j)
		{
			Flag[i][j]	= new vector<size_t> [Nz+1];
			Flagt[i][j]	= new vector<size_t> [Nz+1];

			for (int k=0; k<=Nz; ++k)
			{
				Flag[i][j][k].resize(0);
				Flagt[i][j][k].resize(0);
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

		Flagt[0 ][j][k].push_back(0);
		Flagt[Nx][j][k].push_back(1);
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

		Flagt[i][0 ][k].push_back(2);
		Flagt[i][Ny][k].push_back(3);
	}		
	Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
	Lp[2]->ID = 2;
	Lp[2]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
	Lp[3]->ID = 3;
	Lp[3]->Type = -1;
	if (D==3)
	{
		// For z axis
		for (int i=0; i<=Nx; ++i)
		for (int j=0; j<=Ny; ++j)
		{
			Flag[i][j][0 ].push_back(4);
			Flag[i][j][Nz].push_back(5);

			Flagt[i][j][0 ].push_back(4);
			Flagt[i][j][Nz].push_back(5);
		}
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
	// UpdateFlag(Lp[Lp.size()-1]);
}

inline void DEM::AddDisk2D(int tag, double r, Vector3d& x, double rho)
{
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
	Lp[Lp.size()-1]->SetDisk2D(r);
	// UpdateFlag(Lp[Lp.size()-1]);
}

// inline void DEM::Add2DPolynomialParticle(int tag, VectorXd& coef, Vector3d& x, double rho)
// {
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	Lp[Lp.size()-1]->ID = Lp.size()-1;
// 	Lp[Lp.size()-1]->Set2DPolynomialParticle(coef);
// 	Np += Lp[Lp.size()-1]->P0.size();
// 	Nf += Lp[Lp.size()-1]->Faces.size();
// 	// UpdateFlag(Lp[Lp.size()-1]);
// }

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

inline void DEM::AddTetrahedron(int tag, vector<Vector3d> ver, double rho)
{
	if (ver.size()!=4)
	{
		cout << "\033[1;31mError: Vertices number is not 4, you sure it's tetrahedron?\033[0m\n";		
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
	Vector3d x (0,0,0);
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
	cout << "SetTetrahedron" << endl;
	Lp[Lp.size()-1]->SetTetrahedron(ver);
	cout << "SetTetrahedron Finish" << endl;
}

// inline void DEM::AddDisk2D(int tag, double r, Vector3d& x, double rho)
// {
// 	if (x(2)!=0.)
// 	{
// 		cout << "\033[1;31mError: Disk only works on 2D, check position again.\033[0m\n";		
// 		exit(0);		
// 	}

// 	for (size_t i=0; i<Lp.size(); ++i)
// 	{
// 		if (tag==Lp[i]->Tag)
// 		{
// 			cout << "\033[1;31mError: Cannot add sphere. The tag is existed.\033[0m\n";		
// 			exit(0);
// 		}
// 	}
// 	Lp.push_back(new DISK2D(tag, r, x, rho));
// 	Lp[Lp.size()-1]->ID = Lp.size()-1;
// }

inline void DEM::AddNSpheres(int tag, int np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho)
{
    srand (time(NULL));
    int count = 0;

    Vector3d l = x1-x0;
    cout << l.transpose() << endl;
    for (int n=0; n<1.0e200; n++)
    {
        Vector3d pos(0.,0.,0.);
        for (size_t d=0; d<D; ++d)
        {
        	pos(d) = (l(d)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        }
        // pos(0) = (l(0)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        // pos(1) = (l(1)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        // pos(2) = (l(2)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
        pos += x0;
        bool boolContact = false;

        if (Lp.size()>6)
        {
            // #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t p=6; p<Lp.size(); ++p)
            {
                Vector3d xpq = Lp[p]->X - pos;
                double pqDis = xpq.norm() - r - Lp[p]->R;
                // cout << "pqDis: " << pqDis << endl;
                // cout << "Lp[p]->X: " << Lp[p]->X.transpose() << endl;
                // cout << "pos: " << pos.transpose() << endl;
                if (pqDis<surDis)
                {
                	boolContact = true;
                	break;
                }
            }
        }
        if (boolContact==false)   
        {
            if (D==3)		AddSphere(tag, r, pos, rho);
            else if (D==2)	AddDisk2D(tag, r, pos, rho);
            count++;
            cout << "Already add " << count << " particles" << endl;
            // for (size_t p=2*D; p<Lp.size()-1; ++p)
            // {
            //     Vector3d xpq = Lp[Lp.size()-1]->X - Lp[p]->X;
            //     double pqDis = xpq.norm() - Lp[Lp.size()-1]->R - Lp[p]->R;
            //     if (pqDis<surDis)
            //     {
            //     	cout << "contacted" << endl;
            //     	cout << pqDis << endl;
            //     	abort();
            //     }           
            // }
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
		if (p0->Group==-1)
		{
			if 		(p0->X(0)>Nx)	p0->X(0) = p0->X(0)-Nx-1;
			else if (p0->X(0)<0.)	p0->X(0) = p0->X(0)+Nx+1;
			if 		(p0->X(1)>Ny)	p0->X(1) = p0->X(1)-Ny-1;
			else if (p0->X(1)<0.)	p0->X(1) = p0->X(1)+Ny+1;
			if 		(p0->X(2)>Nz)	p0->X(2) = p0->X(2)-Nz-1;
			else if (p0->X(2)<0.)	p0->X(2) = p0->X(2)+Nz+1;
			p0->VelocityVerlet(Dt);
			p0->UpdateBox(D);
		}
		else
		{
			// #pragma omp critical
			// {
			// 	Lg[p0->Group]->Fh += p0->Fh;
			// 	Lg[p0->Group]->Fc += p0->Fc;
			// }
            for (size_t d=0; d<D; ++d)
            {
                #pragma omp atomic
                Lg[p0->Group]->Fh(d) += p0->Fh(d);
                #pragma omp atomic
                Lg[p0->Group]->Fc(d) += p0->Fc(d);
            }
		}
	}
	// cout << "Lg.size()= " << Lg.size() << endl;
	// For groups, WORNING: rotation is not considered yet.
	for (size_t i=0; i<Lg.size(); ++i)
	{
		DEM_PARTICLE* g0 = Lg[i];
		if 		(g0->X(0)>Nx)	g0->X(0) = g0->X(0)-Nx-1;
		else if (g0->X(0)<0.)	g0->X(0) = g0->X(0)+Nx+1;
		if 		(g0->X(1)>Ny)	g0->X(1) = g0->X(1)-Ny-1;
		else if (g0->X(1)<0.)	g0->X(1) = g0->X(1)+Ny+1;
		if 		(g0->X(2)>Nz)	g0->X(2) = g0->X(2)-Nz-1;
		else if (g0->X(2)<0.)	g0->X(2) = g0->X(2)+Nz+1;
		g0->VelocityVerlet(Dt);
		g0->UpdateBox(D);
		Vector3d displace = g0->X - g0->Xb;
		// cout << "g0->V= " << g0->V.transpose() << endl;

		for (size_t l=0; l<g0->Lp.size(); ++l)
		{
			size_t p = g0->Lp[l];
			DEM_PARTICLE* p0 = Lp[p];
			p0->X += displace;
			p0->V = g0->V;
			if 		(p0->X(0)>Nx)	p0->X(0) = p0->X(0)-Nx-1;
			else if (p0->X(0)<0.)	p0->X(0) = p0->X(0)+Nx+1;
			if 		(p0->X(1)>Ny)	p0->X(1) = p0->X(1)-Ny-1;
			else if (p0->X(1)<0.)	p0->X(1) = p0->X(1)+Ny+1;
			if 		(p0->X(2)>Nz)	p0->X(2) = p0->X(2)-Nz-1;
			else if (p0->X(2)<0.)	p0->X(2) = p0->X(2)+Nz+1;
			p0->UpdateBox(D);
		}
	}
}

inline void DEM::ZeroForceTorque(bool h, bool c)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		Lp[i]->ZeroForceTorque(h, c);
		if ((i-6)<Lg.size())	Lg[i-6]->ZeroForceTorque(h, c);
	}
}

inline void DEM::SetG(Vector3d& g)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		Lp[i]->SetG(g);
		if ((i-6)<Lg.size())	Lg[i-6]->SetG(g);
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

void DEM::LinearDampingPara0(double& kn, double& me, double& gn, double& gt)
{
	// "Introduction to Discrete Element Methods" (Luding)
	// eq(4-5)
	gn 	= -2.*Beta*sqrt(kn*me);
	gt 	= RatioGnt*gn;
}

void DEM::HertzDampingPara0(double& kn, double& me, double& gn, double& gt)
{
	// eq(7) "Granular flow down an inclined plane: Bagnold scaling and rheology"
	// "Micromechanical Origin of Particle Size Segregation"
	gn = -log(Cr)*me;
	gt = RatioGnt*gn;
}

void DEM::DampingParaDoNothing(double& kn, double& me, double& gn, double& gt)
{
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

	(this->*DampingPara)(kn, me, gn, gt);
}

inline void DEM::SetLubrication(double hn, double viscosity)
{
	Hn = hn;
	Viscosity = viscosity;
}

// Contact force model for spheres
inline void DEM::Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, Vector3d& xir, bool& contacted)
{
	// cout << "contact start" << endl;
	contacted = false;
	Vector3d Xi = pi->X;
	Vector3d Xj = pj->X;

	if (pi->ID<6)
	{
		Xi = Xj;
		int axis = pi->ID/2;
		int dirc = pi->ID-2*axis;
		Xi(axis) = dirc*DomSize[axis];
	}

	Vector3d n = Xi-Xj;								// Normal direction (pj pinnts to pi)
	Vector3d cpi, cpj;								// closest point on i and j, contact point
	cpi = Xi; cpj = Xj;								// set to sphere center for sphere collisions
	if (pi->Type==3 && pj->Type==3)					// For polyhedron collisions
	{
		FindClosestPoints3D(pi->P, pj->P, n, cpi, cpj);	// Find closest points
		n = cpi-cpj;									// Contact normal
		// cout << "cpi: " << cpi.transpose() << endl;
		// cout << "cpj: " << cpj.transpose() << endl;
	}
	double delta = pi->R+pj->R-n.norm(); 			// Overlapping distance
	n.normalize();									// Normalize contact normal
	Vector3d cp = cpj+(pj->R-0.5*delta)*n;			// Contact point

	if (pi->crossingFlag || pj->crossingFlag)
	{
		for (size_t i=0; i<pi->Xmir.size(); ++i)
		for (size_t j=0; j<pj->Xmir.size(); ++j)
		{
			double deltat = pi->R+pj->R-(pi->Xmir[i]-pj->Xmir[j]).norm();
			if (deltat>delta)
			{
				delta = deltat;
				n = pi->Xmir[i]-pj->Xmir[j];
			}
		}
	}
	// cout << "delta: " << delta << endl;
	// abort();
	// report contact for updating friction map 
	if (delta>0.)	contacted = true;

	if (contacted)
	{
		// if (delta>0.01*(pi->R+pj->R) && (pi->Group==-1 && pj->Group==-1))
		// {
		// 	cout << "max delta= " << 0.01*(pi->R+pj->R) << endl;
		// 	cout << "xi= " << Xi << endl;
		// 	cout << "xj= " << Xj << endl;
		// 	cout << "delta= " << delta << endl;
		// 	cout << pi->ID << endl;
		// 	cout << pj->ID << endl;
		// 	cout << Xi.transpose() << endl;
		// 	cout << pi->V.transpose() << endl;
		// 	cout << pj->V.transpose() << endl;
		// 	abort();			
		// }
		double kn, gn, kt, gt;
		(this->*ContactPara)(pi, pj, delta, kn, gn, kt, gt);
		Vector3d vn = (pj->V-pi->V).dot(n)*n;			// Relative velocity in normal direction
		Vector3d fn= kn*delta*n + gn*vn;				// Normal contact force
		Vector3d ft (0., 0., 0.);
		Friction(pi, pj, delta, kt, gt, n, fn, xi, ft);	// Friction force and torque
		// Rolling resistance torque
		// double kr = RatioKnr*kn;
		// double gr = RatioGnr*gn;
		// RollingResistance(pi, pj, delta, kr, gr, n, fn, xir, armr);
		// arm += armr;
		Vector3d fnt = fn+ft;							// Total force
		#pragma omp critical
		{
			pi->Fc += fnt;
			pj->Fc -= fnt;

			pi->Tc += (pi->Qfi._transformVector(fnt)).cross(pi->Qfi._transformVector(Xi-cp));
			pj->Tc -= (pj->Qfi._transformVector(fnt)).cross(pj->Qfi._transformVector(Xj-cp));
		}
   //      for (size_t d=0; d<D; ++d)
   //      {
   //          #pragma omp atomic
   //      	pi->Fc(d) += fn(d)+ft(d);
   //          #pragma omp atomic
   //          pj->Fc(d) -= fn(d)+ft(d);
   //          #pragma omp atomic
			// pi->Tc(d) += (pi->R-0.5*delta)*arm(d);
			// #pragma omp atomic
			// pj->Tc(d) += (pj->R-0.5*delta)*arm(d);
   //      }		
	}
}

inline void DEM::Friction(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kt, double gt, Vector3d& n, Vector3d& fn, Vector3d& xi, Vector3d& ft)
{
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
	// Static tangential force
	double fts = FsTable[pi->MID][pj->MID]*fn.norm();
	// Tangential force
	ft = -kt*xi-gt*vt;									// eq.18
	if (ft.norm()>fts)
	{
		Vector3d t = ft.normalized();
		ft = fts*t;										// eq.21
		xi = -(fts*t+gt*vt)/kt;							// eq.20
	}
}

inline void DEM::RollingResistance(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kr, double gr, Vector3d& n, Vector3d& fn, Vector3d& xir, Vector3d& armr)
{
	// Reduced radius
	double rij = pi->R*pj->R/(pi->R+pj->R);
	// Rolling velocity
	Vector3d vij = rij*(n.cross(pj->W)-n.cross(pi->W));	// eq.15
	// Relative tangential velocity at the contact point
	Vector3d vt = vij-n.dot(vij)*n;						// eq.9
	// Update tangential spring
	Vector3d xi0 = vt*Dt;
	size_t key = Key(pi->ID, pj->ID);
	if (RMap.count(key)>0)	xi0 += RMap[key];			// eq.19
	// Project to current tangential plane
	xir = xi0 - n.dot(xi0)*n;							// eq.17
	// Static tangential force
	double fts = RsTable[pi->MID][pj->MID]*fn.norm();
	// Tangential force
	Vector3d ft = -kr*xir-gr*vt;						// eq.18
	if (ft.norm()>fts)
	{
		Vector3d t = ft.normalized();
		ft = fts*t;										// eq.21
		xir = -(fts*t+gr*vt)/kr;						// eq.20
	}
	// Torque with normalized arm
	armr = -n.cross(ft);
}

inline double DEM::EffectiveValue(double ai, double aj)
{
	double a = 0.;
	if (ai>1.0e-12 && aj>1.0e-12)	a = ai*aj/(ai+aj);
	return (a);
}

// Update Xmir for periodic BC
inline void DEM::UpdateXmir(DEM_PARTICLE* p0)
{
	// cout << p0->X.transpose() << endl;
	p0->Xmir.clear();
	// add original position
	Vector3d xmir0 = p0->X;
	p0->Xmir.push_back(xmir0);
	// cout << "1111111111111111" << endl;
	// check crossing first
	if (p0->crossingFlag)
	{
		// cout << "crossingFlag " << endl;
		// cout << p0->X.transpose() << endl;
		// list of periodic dircetion
		vector<size_t> dirc;
		for (size_t i=0; i<D; ++i)
		{
			if (p0->crossing[i])	dirc.push_back(i);
		}
		// add mirror position for first direction
		size_t d0 = dirc[0];
		Vector3d xmir = p0->X;
		if (xmir(d0)+p0->R>DomSize[d0])		xmir(d0) -= DomSize[d0]+1;
		else if (xmir(d0)-p0->R<0.)			xmir(d0) += DomSize[d0]+1;
		p0->Xmir.push_back(xmir);
		// cout << "222222222222222" << endl;
		// add mirror position for second direction
		if (dirc.size()>1)
		{
			size_t d1 = dirc[1];
			// for every exiting mirror positions add a new one for this direction 2->4
			for (size_t i=0; i<2; ++i)
			{
				Vector3d xmiri = p0->Xmir[i];
				if (xmiri(d1)+p0->R>DomSize[d1])	xmiri(d1) -= DomSize[d1]+1;
				else if (xmiri(d1)-p0->R<0.)		xmiri(d1) += DomSize[d1]+1;
				p0->Xmir.push_back(xmiri);
			}
		}
		// cout << "3333333333333" << endl;
		// add mirror position for third direction
		if (dirc.size()>2)
		{
			size_t d2 = dirc[2];
			// for every exiting mirror positions add a new one for this direction 4->8
			for (size_t i=0; i<4; ++i)
			{
				Vector3d xmiri = p0->Xmir[i];
				if (xmiri(d2)+p0->R>DomSize[d2])	xmiri(d2) -= DomSize[d2]+1;
				else if (xmiri(d2)-p0->R<0.)		xmiri(d2) += DomSize[d2]+1;
				p0->Xmir.push_back(xmiri);
			}
		}
		// cout << "4444444444444444" << endl;
	}
}

inline void DEM::UpdateXmirGlobal()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		UpdateXmir(Lp[i]);
	}	
}

inline void DEM::UpdateFlag(DEM_PARTICLE* p0)
{
	for (int i=p0->Min(0); i<=p0->Max(0); ++i)
	for (int j=p0->Min(1); j<=p0->Max(1); ++j)
	for (int k=p0->Min(2); k<=p0->Max(2); ++k)
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

// inline void DEM::FindContactBasedOnNode(bool first)
// {
// 	// cout << "init flag" << endl;
// 	if (first)
// 	{
// 		for (size_t p=6; p<Lp.size(); ++p)
// 		{
// 			DEM_PARTICLE* p0 = Lp[p];
// 			/*if (!p0->fixed)	*/UpdateFlag(p0);
// 		}		
// 	}
// 	CMap.clear();
// 	// cout << "update flagt" << endl;
// 	// update Flagt based on its neighbours
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (int i=0; i<=Nx; ++i)
//     for (int j=0; j<=Ny; ++j)
//     for (int k=0; k<=Nz; ++k)
//     {
//     	Vector3i ind (i,j,k);
//     	Vector3d indd (i,j,k);
//     	// cout << "ind: " << ind.transpose() << endl;
//     	// a list of potental particle ID which may mark this cell 
//     	vector<size_t> li;
//     	for (size_t n=0; n<Ne.size(); ++n)
//     	{
//     		Vector3i c = ind+Ne[n];
//     		if (c(0)>-1 && c(0)<=Nx && c(1)>-1 && c(1)<=Ny && c(2)>-1 && c(2)<=Nz)
//     		{
// 				li.insert(li.end(), Flag[c(0)][c(1)][c(2)].begin(), Flag[c(0)][c(1)][c(2)].end());
//     		}
//     	}
//     	sort( li.begin(), li.end() );
//     	li.erase( unique( li.begin(), li.end() ), li.end() );

//     	for (size_t l=0; l<li.size(); ++l)
//     	{
// 			// Distance from the cell centre to particle surface
// 			double dis = 0.;
// 			size_t p = li[l];
// 			if (Lp[p]->Type==1)	dis = (indd-Lp[p]->X).norm()-Lp[p]->R-0.87;
// 			if (dis<0)
// 			{
// 				Flagt[i][j][k].push_back(p);
// 			}
//     	}
//     	// if at least 2 particles marked this cell
//     	if (Flagt[i][j][k].size()>1)
//     	{
//     		for (size_t m=0; m<Flagt[i][j][k].size(); ++m)
//     		for (size_t o=m+1; o<Flagt[i][j][k].size(); ++o)
//     		{
// 				size_t p = Flagt[i][j][k][m];
// 				size_t q = Flagt[i][j][k][o];

// 				size_t min0 = min(p,q);
// 				size_t max0 = max(p,q);
// 				size_t key = Key(min0,max0);
// 				#pragma omp critical
// 				{
// 					if (!CMap[key])
// 					{
// 						Lc.push_back({min0, max0});
// 						CMap[key] = true;
// 					}
// 				}
//     		}
//     	}
//     }
//     // cout << "reset flag" << endl;
//     // reset Flag to Flagt
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (int i=0; i<=Nx; ++i)
//     for (int j=0; j<=Ny; ++j)
//     for (int k=0; k<=Nz; ++k)
//     {
//     	bool show = false;
//     	Flag[i][j][k] = Flagt[i][j][k];
//     	// if (Flag[i][j][k].size()>0)	show = true;
//     	// if (show)	cout << Flag[i][j][k].size() << endl;
//     	Flagt[i][j][k].clear();
//     	// if (show)	cout << Flag[i][j][k].size() << endl;
//     	// if (show)	cout << Flagt[i][j][k].size() << endl;
//     	// if (show) abort();
//     }
//     // cout << "done flag" << endl;
// }

inline void DEM::FindContact()
{
	CMap.clear();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];

		for (int i=p0->Min(0); i<=p0->Max(0); ++i)
		for (int j=p0->Min(1); j<=p0->Max(1); ++j)
		for (int k=p0->Min(2); k<=p0->Max(2); ++k)
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
					if (/*Lp[ind]->fixed ||*/ ind<6)
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
		/*if (!p0->fixed)	*/UpdateFlag(p0);
	}
}

inline void DEM::Contact(bool writeFc, int n)
{
    if (Lc.size()>0.)
    {
    	unordered_map<size_t, Vector3d> fmap;
    	unordered_map<size_t, Vector3d> rmap;
    	// #pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t l=0; l<Lc.size(); ++l)
		{
			int i = Lc[l][0];
			int j = Lc[l][1];
			// if (i==2)
			// {
			// 	cout << "i= " << i << " j= " << j << endl;
			// 	// abort();
			// }
			// cout << "i= " << i << " j= " << j << endl;
			// abort();
			bool contacted = false;
			Vector3d xi (0.,0.,0.);
			Vector3d xir (0.,0.,0.);
			Contact2P(Lp[i], Lp[j], xi, xir, contacted);
			if (contacted)
			{
				fmap[Key(i,j)] = xi;
				rmap[Key(i,j)] = xir;
			}
		}
		FMap = fmap;
		RMap = rmap;
    }
    if (writeFc)	WriteContactForceFileH5(n);
	// Lc.clear();
}

inline void DEM::Solve(int tt, int ts, double dt, bool writefile)
{
	Dt = dt;
	for (int t=0; t<tt; ++t)
	{
		bool show = false;
		if (t%100==0)	show = true;
		if (t%ts == 0)
		{
			cout << "Time Step ============ " << t << endl;
			if (writefile /*&& t>4000*/)	WriteFileH5(t);
		}

		if (show)	cout << "Time Step ============ " << t << endl;
		clock_t t_start = std::clock();
		FindContact();
		Lc.push_back({6,7});
		// bool firstStep = false;
		// if (t==0)	firstStep = true;
		// FindContactBasedOnNode(firstStep);
		clock_t t_end = std::clock();
		if (show)	cout << "FindContact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::clock();
		Contact(false, 0);
		int nc = Lc.size();
		Lc.clear();
		t_end = std::clock();
		if (show)   cout << "nc= " << nc << endl;
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

		// double ta = 1e5;
		// // Vector3d vtop (0.008*sin(t/ta), 0., 0.);
		// // if (t>ta)	vtop(0)=0.008;
		// double v = 0.003;
		// Vector3d vtop (v*(t/ta), 0., 0.);
		// if (vtop(0)>v)	vtop(0)=v;
		// Vector3d wtop (0., 0., 0.);
		// for (size_t i=0; i<Lp.size(); ++i)
		// {
		// 	if (Lp[i]->Tag==1)
		// 	{
		// 		Lp[i]->FixV(vtop);
		// 		Lp[i]->FixW(wtop);
		// 	}
		// }

		// for (size_t p=7; p<Lp.size(); ++p)
		// {
		// 	DEM_PARTICLE* p0 = Lp[p];

		// 	double delta = p0->R+Lp[6]->R-(p0->X-Lp[6]->X).norm();
		// 	size_t key = Key(6,p);
		// 	if (delta>1.0e-4 && ! CMap[key])
		// 	{
		// 		cout << "p= " << p << endl;
		// 		// size_t key = Key(6,p);
		// 		cout << "CMap[key]= " << CMap[key] << endl;
		// 		cout << "p0->X= " << p0->X << endl;
		// 		cout << "Lp[6]->X= " << Lp[6]->X << endl;
		// 		abort();
		// 	}
		// }

		Move();
		t_end = std::clock();
		if (show)	cout << "Move time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::clock();
		// WriteFileParticleInfo(t);
		ZeroForceTorque(true, true);
		t_end = std::clock();
		if (show)	cout << "ZeroForceTorque time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	}
}

void DEM::DeleteParticles()
{
	vector <DEM_PARTICLE*>	Lpt;
	Lpt.resize(0);

	for (size_t p=0; p<Lp.size(); ++p)
	{
		if (!Lp[p]->removed)	Lpt.push_back(Lp[p]);
	}
	Lp = Lpt;

	for (size_t p=0; p<Lp.size(); ++p)
	{
		Lp[p]->ID = p;
	}
}

void DEM::LoadDEMFromH5( string fname, double scale, double rhos)
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

	// H5std_string DATASET_NAME_RHO( "Rho" );
	// H5File file_rho( FILE_NAME, H5F_ACC_RDONLY );
	// DataSet dataset_rho = file_rho.openDataSet( DATASET_NAME_RHO );
	// DataSpace dataspace_rho = dataset_r.getSpace();
 //    hsize_t dims_rho[2];
 //    dataspace_rho.getSimpleExtentDims( dims_rho, NULL);
 //    hsize_t dimsm_rho = dims_rho[0];

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

    // double data_rho[dimsm_rho];
    // dataset_rho.read( data_rho, PredType::NATIVE_DOUBLE, dataspace_rho, dataspace_rho );

    double data_tag[dimsm_tag];
    dataset_tag.read( data_tag, PredType::NATIVE_DOUBLE, dataspace_tag, dataspace_tag );

    int np = dimsm_pos/3;
    for (int i=6; i<np; ++i)
    {
    	Vector3d pos (scale*data_pos[3*i], scale*data_pos[3*i+1], scale*data_pos[3*i+2]);
    	double r = scale*data_r[i];
    	int tag = (int) data_tag[i];
    	// double rho = data_rho[i];
    	double rho = rhos;
    	AddSphere(tag, r, pos, rho);
    }
    cout << "========= Loaded "<< Lp.size()-6<< " DEM particles from " << fname << "==============" << endl;
}

// inline void DEM::WriteFileH5(int n)
// {
// 	stringstream	out;							//convert int to string for file name.
// 	out << setw(9) << setfill('0') << n;			
// 	string file_name_h5 = "DEM_"+out.str()+".h5";

// 	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
// 	hsize_t	dims_scalar[1] = {Lp.size()};			//create data space.
// 	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.

// 	hsize_t n_points = 8*Lp.size();
// 	hsize_t n_faces  = 6*Lp.size();
// 	hsize_t	dims_points [1] = {3*n_points};			//create data space.
// 	hsize_t	dims_faces  [1] = {5*n_faces};			//create data space.
// 	hsize_t	dims_fscalar[1] = {n_faces};			//create data space.

// 	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
// 	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

// 	int rank_points = sizeof(dims_points ) / sizeof(hsize_t);
// 	int rank_faces  = sizeof(dims_faces  ) / sizeof(hsize_t);
// 	int rank_fscalar= sizeof(dims_fscalar) / sizeof(hsize_t);

// 	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
// 	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

// 	DataSpace	*space_points = new DataSpace(rank_points, dims_points);
// 	DataSpace	*space_faces  = new DataSpace(rank_faces , dims_faces );
// 	DataSpace	*space_fscalar= new DataSpace(rank_fscalar,dims_fscalar);

// 	double* r_h5 	= new double[  Lp.size()];
// 	double* rho_h5 	= new double[  Lp.size()];
// 	double* tag_h5 	= new double[  Lp.size()];
// 	double* pos_h5 	= new double[3*Lp.size()];
// 	double* vel_h5 	= new double[3*Lp.size()];
// 	double* agv_h5	= new double[3*Lp.size()];
// 	double* fh_h5 	= new double[3*Lp.size()];

// 	double* poi_h5	= new double[3*n_points];
// 	int* 	fac_h5	= new int   [5*n_faces];

// 	double* fv_h5 	= new double[n_faces];

// 	int count_p = 0;
// 	int count_f = 0;
// 	int count_fv = 0;

// 	for (size_t i=0; i<Lp.size(); ++i)
// 	{
// 		Vector3d agv = Lp[i]->Q._transformVector(Lp[i]->W);
//         r_h5  [  i  ] 	= Lp[i]->R;
//         rho_h5[  i  ] 	= Lp[i]->Rho;
//         tag_h5[  i  ] 	= Lp[i]->Tag;
// 		pos_h5[3*i  ] 	= Lp[i]->X(0);
// 		pos_h5[3*i+1] 	= Lp[i]->X(1);
// 		pos_h5[3*i+2] 	= Lp[i]->X(2);
// 		vel_h5[3*i  ] 	= Lp[i]->V(0);
// 		vel_h5[3*i+1] 	= Lp[i]->V(1);
// 		vel_h5[3*i+2] 	= Lp[i]->V(2);
// 		fh_h5[3*i  ] 	= Lp[i]->Fh(0);
// 		fh_h5[3*i+1] 	= Lp[i]->Fh(1);
// 		fh_h5[3*i+2] 	= Lp[i]->Fh(2);
// 		agv_h5[3*i  ] 	= agv(0);
// 		agv_h5[3*i+1] 	= agv(1);
// 		agv_h5[3*i+2] 	= agv(2);

// 		for (size_t j=0; j<Lp[i]->P.size(); ++j)
// 		{
// 			poi_h5[count_p  ] = Lp[i]->P[j](0);
// 			poi_h5[count_p+1] = Lp[i]->P[j](1);
// 			poi_h5[count_p+2] = Lp[i]->P[j](2);
// 			count_p += 3;
// 		}
// 		// cout << "Lp[i]->Faces.size()= " << Lp[i]->Faces.size() << endl;
// 		for (size_t k=0; k<Lp[i]->Faces.size(); ++k)
// 		{
// 			fac_h5[count_f  ] = 5;
// 			fac_h5[count_f+1] = Lp[i]->Faces[k](0);
// 			fac_h5[count_f+2] = Lp[i]->Faces[k](1);
// 			fac_h5[count_f+3] = Lp[i]->Faces[k](2);
// 			fac_h5[count_f+4] = Lp[i]->Faces[k](3);
// 			count_f += 5;
// 			fv_h5[count_fv] = Lp[i]->V.norm();
// 			count_fv++;
// 		}
// 	}

// 	DataSet	*dataset_r 		= new DataSet(file.createDataSet("Radius", PredType::NATIVE_DOUBLE, *space_scalar));
// 	DataSet	*dataset_rho	= new DataSet(file.createDataSet("Rho", PredType::NATIVE_DOUBLE, *space_scalar));
// 	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
//     DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
//     DataSet	*dataset_agv	= new DataSet(file.createDataSet("AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector));
//     DataSet	*dataset_fh		= new DataSet(file.createDataSet("HydroForce", PredType::NATIVE_DOUBLE, *space_vector));

//     DataSet	*dataset_poi	= new DataSet(file.createDataSet("Points", PredType::NATIVE_DOUBLE, *space_points));
//     DataSet	*dataset_fac	= new DataSet(file.createDataSet("Faces", PredType::NATIVE_INT, *space_faces));
//     DataSet	*dataset_fv		= new DataSet(file.createDataSet("FaceVelocity", PredType::NATIVE_DOUBLE, *space_fscalar));

// 	dataset_r->write(r_h5, PredType::NATIVE_DOUBLE);
// 	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
// 	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
// 	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
// 	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
// 	dataset_agv->write(agv_h5, PredType::NATIVE_DOUBLE);
// 	dataset_fh->write(fh_h5, PredType::NATIVE_DOUBLE);

// 	dataset_poi->write(poi_h5, PredType::NATIVE_DOUBLE);
// 	dataset_fac->write(fac_h5, PredType::NATIVE_INT);
// 	dataset_fv->write(fv_h5, PredType::NATIVE_DOUBLE);

// 	delete space_scalar;
// 	delete space_vector;

// 	delete dataset_r;
// 	delete dataset_rho;
// 	delete dataset_tag;
// 	delete dataset_pos;
// 	delete dataset_vel;
// 	delete dataset_agv;
// 	delete dataset_fh;

// 	delete dataset_poi;
// 	delete dataset_fac;
// 	delete dataset_fv;

// 	delete r_h5;
// 	delete rho_h5;
// 	delete tag_h5;
// 	delete pos_h5;
// 	delete vel_h5;
// 	delete agv_h5;
// 	delete fh_h5;

// 	delete poi_h5;
// 	delete fac_h5;
// 	delete fv_h5;

// 	file.close();

// 	string file_name_xmf = "DEM_"+out.str()+".xmf";

//     std::ofstream oss;
//     oss.open(file_name_xmf);
//     oss << "<?xml version=\"1.0\" ?>\n";
//     oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
//     oss << "<Xdmf Version=\"2.0\">\n";
//     oss << " <Domain>\n";

//     if (n_faces>0)
//     {
// 	    oss << "   <Grid Name=\"DEM_FACES\">\n";
// 	    oss << "     <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << n_faces << "\">\n";
// 	    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << 5*n_faces << "\">\n";
// 	    oss << "        " << file_name_h5 <<":/Faces \n";
// 	    oss << "       </DataItem>\n";
// 	    oss << "     </Topology>\n";
// 	    oss << "     <Geometry GeometryType=\"XYZ\">\n";
// 	    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_points << " 3\" >\n";
// 	    oss << "        " << file_name_h5 <<":/Points \n";
// 	    oss << "       </DataItem>\n";
// 	    oss << "     </Geometry>\n";
// 	    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
// 	    oss << "       <DataItem Dimensions=\"" << n_faces << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
// 	    oss << "        " << file_name_h5 <<":/FaceVelocity\n";
// 	    oss << "       </DataItem>\n";
// 	    oss << "     </Attribute>\n";
// 	    oss << "   </Grid>\n";
//     }

//     oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
//     oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Lp.size() << "\"/>\n";
//     oss << "     <Geometry GeometryType=\"XYZ\">\n";
//     oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Lp.size() << " 3\" >\n";
//     oss << "        " << file_name_h5 <<":/Position \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Geometry>\n";
//     oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Radius \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Rho \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Tag \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Velocity\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/AngularVelocity\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"HydroForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/HydroForce\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "   </Grid>\n";
//     oss << " </Domain>\n";
//     oss << "</Xdmf>\n";
//     oss.close();
// }

inline void DEM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(9) << setfill('0') << n;			
	string file_name_h5 = "DEM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	size_t npar = Lp.size()-6;

	hsize_t	dims_scalar[1] = {npar};			//create data space.
	hsize_t	dims_vector[1] = {3*npar};			//create data space.

	hsize_t n_points = 0;
	hsize_t n_faces  = 0;
	hsize_t n_fe  = 0;

	for (size_t i=6; i<Lp.size(); ++i)
	{
		n_points += Lp[i]->P.size();
		n_faces  += Lp[i]->Faces.size();
		n_fe 	 += Lp[i]->Nfe;
	}

	hsize_t	dims_points [1] = {3*n_points};			//create data space.
	hsize_t	dims_faces  [1] = {n_fe};				//create data space.
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

	double* r_h5 	= new double[  npar];
	double* rho_h5 	= new double[  npar];
	double* tag_h5 	= new double[  npar];
	double* pos_h5 	= new double[3*npar];
	double* vel_h5 	= new double[3*npar];
	double* agv_h5	= new double[3*npar];
	double* fh_h5 	= new double[3*npar];

	double* poi_h5	= new double[3*n_points];
	int* 	fac_h5	= new int   [n_fe];

	double* fv_h5 	= new double[n_faces];
	double* ftag_h5 = new double[n_faces];

	int count_p = 0;
	int count_f = 0;
	int count_fv = 0;

	for (size_t i=0; i<npar; ++i)
	{
		size_t ind = i+6;
		Vector3d agv = Lp[ind]->Q._transformVector(Lp[ind]->W);
        r_h5  [  i  ] 	= Lp[ind]->R;
        rho_h5[  i  ] 	= Lp[ind]->Rho;
        tag_h5[  i  ] 	= Lp[ind]->Tag;
		pos_h5[3*i  ] 	= Lp[ind]->X(0);
		pos_h5[3*i+1] 	= Lp[ind]->X(1);
		pos_h5[3*i+2] 	= Lp[ind]->X(2);
		vel_h5[3*i  ] 	= Lp[ind]->V(0);
		vel_h5[3*i+1] 	= Lp[ind]->V(1);
		vel_h5[3*i+2] 	= Lp[ind]->V(2);
		fh_h5[3*i  ] 	= Lp[ind]->Fh(0);
		fh_h5[3*i+1] 	= Lp[ind]->Fh(1);
		fh_h5[3*i+2] 	= Lp[ind]->Fh(2);
		agv_h5[3*i  ] 	= agv(0);
		agv_h5[3*i+1] 	= agv(1);
		agv_h5[3*i+2] 	= agv(2);

		for (size_t j=0; j<Lp[ind]->P.size(); ++j)
		{
			poi_h5[count_p  ] = Lp[ind]->P[j](0);
			poi_h5[count_p+1] = Lp[ind]->P[j](1);
			poi_h5[count_p+2] = Lp[ind]->P[j](2);
			count_p += 3;
		}
		size_t add_poi = i*Lp[ind-1]->P.size();
		for (size_t k=0; k<Lp[ind]->Faces.size(); ++k)
		{
			size_t num = Lp[ind]->Faces[k].size()+1;
			fac_h5[count_f  ] = num;
			for (size_t m=0; m<num-1; ++m)
			{
				fac_h5[count_f+m+1] = Lp[ind]->Faces[k](m)+add_poi;
			}

			count_f += num;
			fv_h5[count_fv] 	= Lp[ind]->V.norm();
			ftag_h5[count_fv] 	= Lp[ind]->Tag;
			count_fv++;
		}
	}

	DataSet	*dataset_r 		= new DataSet(file.createDataSet("Radius", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_rho	= new DataSet(file.createDataSet("Rho", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_agv	= new DataSet(file.createDataSet("AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_fh		= new DataSet(file.createDataSet("HydroForce", PredType::NATIVE_DOUBLE, *space_vector));

    DataSet	*dataset_poi	= new DataSet(file.createDataSet("Points", PredType::NATIVE_DOUBLE, *space_points));
    DataSet	*dataset_fac	= new DataSet(file.createDataSet("Faces", PredType::NATIVE_INT, *space_faces));
    DataSet	*dataset_fv		= new DataSet(file.createDataSet("FaceVelocity", PredType::NATIVE_DOUBLE, *space_fscalar));
    DataSet	*dataset_ftag	= new DataSet(file.createDataSet("FaceTag", PredType::NATIVE_DOUBLE, *space_fscalar));

	dataset_r->write(r_h5, PredType::NATIVE_DOUBLE);
	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
	dataset_agv->write(agv_h5, PredType::NATIVE_DOUBLE);
	dataset_fh->write(fh_h5, PredType::NATIVE_DOUBLE);

	dataset_poi->write(poi_h5, PredType::NATIVE_DOUBLE);
	dataset_fac->write(fac_h5, PredType::NATIVE_INT);
	dataset_fv->write(fv_h5, PredType::NATIVE_DOUBLE);
	dataset_ftag->write(ftag_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;
	delete dataset_r;
	delete dataset_rho;
	delete dataset_tag;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_agv;
	delete dataset_fh;
	delete dataset_poi;
	delete dataset_fac;
	delete dataset_fv;
	delete dataset_ftag;
	delete r_h5;
	delete rho_h5;
	delete tag_h5;
	delete pos_h5;
	delete vel_h5;
	delete agv_h5;
	delete fh_h5;
	delete poi_h5;
	delete fac_h5;
	delete fv_h5;
	delete ftag_h5;

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
	    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << n_fe << "\">\n";
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
	    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
	    oss << "       <DataItem Dimensions=\"" << n_faces << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
	    oss << "        " << file_name_h5 <<":/FaceTag\n";
	    oss << "       </DataItem>\n";
	    oss << "     </Attribute>\n";
	    oss << "   </Grid>\n";
    }

    oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << npar << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << npar << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Rho \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/AngularVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"HydroForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/HydroForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}
inline void DEM::WriteContactForceFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(9) << setfill('0') << n;			
	string file_name_h5 = "DEM_Force_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lc.size()};			//create data space.
	hsize_t	dims_vector[1] = {3*Lc.size()};			//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* p0_h5 	= new double[  Lc.size()];
	double* p1_h5 	= new double[  Lc.size()];
	double* fc_h5 	= new double[3*Lc.size()];

	for (size_t i=0; i<Lc.size(); ++i)
	{
        p0_h5[  i  ] 	= Lc[i][0];
        p1_h5[  i  ] 	= Lc[i][1];
		fc_h5[3*i  ] 	= Lp[Lc[i][1]]->Fc(0);
		fc_h5[3*i+1] 	= Lp[Lc[i][1]]->Fc(1);
		fc_h5[3*i+2] 	= Lp[Lc[i][1]]->Fc(2);
	}

	DataSet	*dataset_p0 = new DataSet(file.createDataSet("P0", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_p1	= new DataSet(file.createDataSet("P1", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_fc	= new DataSet(file.createDataSet("ContactForce", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_p0->write(p0_h5, PredType::NATIVE_DOUBLE);
	dataset_p1->write(p1_h5, PredType::NATIVE_DOUBLE);
	dataset_fc->write(fc_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;

	delete dataset_p0;
	delete dataset_p1;
	delete dataset_fc;

	delete p0_h5;
	delete p1_h5;
	delete fc_h5;

	file.close();
}

inline void DEM::WriteFileParticleInfo(int n)
{
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<Lp.size(); ++p)
    {
    	if (p==1139 || p==2162)
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