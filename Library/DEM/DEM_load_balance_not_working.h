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
// #include <UNITSPHERE.h>
// #include <2D_PDEM_FUNCTIONS.h>

// Processor subdomain
struct SUBDOMAIN
{
	size_t Nx;		// subdomain size
	size_t Ny;
	size_t Nz;
	size_t Ncz;		// used to calculate cell index
	size_t Ncy;
	int Nc;			// total number of cells

	Vector3d Xmin;		// min subdomain position
	Vector3d Xmax;		// max subdomain position
	/*===============================================================*/
	/*                    {this part is particles close to boundary} */ 
	/* Lp[p0, p1, p2......|pn......................................|]*/
	/*===============================================================*/
	vector<pair<size_t,bool>> Lp;	// list of particle belongs to this subdomain
	// vector<bool> Lb;	// list of particle status(true for boundary particles)
};

class DEM
{
public:
	DEM(int nx, int ny, int nz, string cmtype, string dmtype, double cr);
	~DEM();
	void Init(Vector3d& lx);
	void AddPlane(int tag, Vector3d& x, Vector3d& n);
	void AddSphere(int tag, double r, Vector3d& x, double rho);
	void AddCuboid(int tag, double lx, double ly, double lz, Vector3d& x, double rho);
	void AddTetrahedron(int tag, vector<Vector3d> ver, double rho);
	void AddDisk2D(int tag, double r, Vector3d& x, double rho);
	void AddNSpheres(int tag, int np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho);
	void AddNMetaballs(int tag, size_t seed, int np, Vector3d& x0, Vector3d& x1, double rho, double rs, vector<Vector3d>& metaP, VectorXd& metaK);
	void Add2DPolynomialParticle(int tag, VectorXd& coef, Vector3d& x, double rho);
	void AddMetaball(int tag, double rs, vector<Vector3d>& metaP, VectorXd& metaK, double rho);
	void Move(bool movePoints);
	void UpdatePositionGlobal(double dt, bool movePoints, bool resetHydro);
	void UpdateVelocityGlobal(double dt);
	void ZeroForceTorque(bool h, bool c);
	void SetG(Vector3d& g);
	double EffectiveValue(double ai, double aj);											// Calculate effective values for contact force
	void RecordX();																			// Record position at Xb for check refilling LBM nodes
	void Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, Vector3d& xir, Vector3d& finalP, bool& contacted);
	// void Friction(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kt, double gt, Vector3d& n, Vector3d& fn, Vector3d& xi, Vector3d& ft);
	void Friction(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double kt, double gt, Vector3d& cp, Vector3d& n, Vector3d& fn, Vector3d& xi, Vector3d& ft);
	void RollingResistance(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kr, double gr, Vector3d& n, Vector3d& fn, Vector3d& xir, Vector3d& armr);
	void UpdateFlag(DEM_PARTICLE* p0);
	void UpdateXmir(DEM_PARTICLE* p0);
	void UpdateXmirGlobal();
	void SetLinkedCell(Vector3d& lx);
	void GenBinsProc();
	void LinkedCell();
	void FindContact();
	void FindContactBasedOnNode(bool first);
	void Contact(bool writeFc, int n);
	void LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void LinearDampingPara0(double& kn, double& me, double& gn, double& gt);
	void HertzDampingPara0(double& kn, double& me, double& gn, double& gt);
	void DampingParaDoNothing(double& kn, double& me, double& gn, double& gt);
	void SetLubrication(double hn, double viscosity);
	void SolveOneStep(double dt, bool save);
	void Solve(int tt, int ts, double dt, bool writefile);
	void DeleteParticles();
	void LoadDEMFromH5( string fname, double scale, double rhos);
	void WriteFileH5(int n);
	void WriteContactForceFileH5(int n);
	void WriteFileParticleInfo(int n);
	void ActiveMetaball();
	void FindIndex(size_t n, size_t& i, size_t& j, size_t& k);

	void (DEM::*ContactPara)(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double& kn, double& gn, double& kt, double& gt);
	void (DEM::*DampingPara)(double& kn, double& me, double& gn, double& gt);

	bool 							Periodic[3];

	size_t 							Nproc;
    size_t 							D;														// Dimension
    int 							Nx;														// Mesh size for contact detection
    int 							Ny;
    int 							Nz;
    int 							Ncz;
    int 							Ncy;
    int 							Nc;														// Total number of cells for linked cell

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

	Vector3d 						Lx;														// Size of cells for linked cell method

	vector<Vector3d> 				Ne;														// Relative location of neighbor cells
	vector<Vector3i> 				Nei;													// Relative location (int) of neighbor cells

	vector<int> 					NeiV;

	vector<Vector3d> 				UnitSpherePoints;										// Unit sphere mesh points for metaball
	vector<VectorXi> 				UnitSphereFaces;										// Unit sphere mesh faces for metaball

    vector<size_t>***		 		Flag;													// Flag of lattice type
    vector<size_t>***		 		Flagt;													// Flag of lattice type

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < DEM_PARTICLE* >		Lg;														// List of groups
	vector < vector< size_t > >     Lc;                                                 	// List of potential contacted paricles' ID
	
	vector < SUBDOMAIN* >			Lsd;													// List of subdomain

	vector<vector<size_t>> 			Bins;													// List of bins for linked cell method
	vector<vector<size_t>>			BinsProc;												// List of actived bins for each proc

	unordered_map<size_t, bool> 	CMap;													// Contact Map
	unordered_map<size_t, Vector3d> FMap;													// Friction Map
	unordered_map<size_t, Vector3d> RMap;													// Rolling resistance Map

	unordered_map<size_t, Vector3d> MMap;												// metaball init point Map
	vector<unordered_map<size_t, Vector3d>> MMapt;										// metaball init point Map

	double 							FsTable[10][10];										// Friction coefficient (static) table
	double 							FdTable[10][10];										// Friction coefficient (dynamic) table
	double 							RsTable[10][10];										// Rolling friction coefficient (static) table
	double 							RdTable[10][10];										// Rolling Friction coefficient (dynamic) table

	double tmpd;
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
		Nei  = {   { 0, 0, 0},
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
	// MMap.clear();
	// MMapt.resize(Nproc);
	tmpd = -1.e20;
}

void DEM::ActiveMetaball()
{
	// LoadUnitSphereMesh(UnitSpherePoints, UnitSphereFaces);
	ReadOBJ("UnitSphere.obj", UnitSpherePoints, UnitSphereFaces);
	cout << "UnitSpherePoints: " << UnitSpherePoints.size() << endl;
	cout << "Metaball actived." << endl;
}

// Map a pair of integer to one integer key for hashing
// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
inline size_t Key(int i, int j)
{
	return (i+j+1)*(i+j)/2+j;
}

// inline void DEM::Init()
// {
//     cout << "================ Start init. ================" << endl;

//     Dt = 1.;

// 	Flag = new vector<size_t>** [Nx+1];
// 	Flagt = new vector<size_t>** [Nx+1];

// 	for (int i=0; i<=Nx; ++i)
// 	{
// 		Flag[i] = new vector<size_t>* [Ny+1];
// 		Flagt[i] = new vector<size_t>* [Ny+1];

// 		for (int j=0; j<=Ny; ++j)
// 		{
// 			Flag[i][j]	= new vector<size_t> [Nz+1];
// 			Flagt[i][j]	= new vector<size_t> [Nz+1];

// 			for (int k=0; k<=Nz; ++k)
// 			{
// 				Flag[i][j][k].resize(0);
// 				Flagt[i][j][k].resize(0);
// 			}
// 		}
// 	}
// 	// Mark flag and add particles for the surounding box to handle walls and periodic BC
// 	// For x axis
// 	for (int j=0; j<=Ny; ++j)
// 	for (int k=0; k<=Nz; ++k)
// 	{
// 		Flag[0 ][j][k].push_back(0);
// 		Flag[Nx][j][k].push_back(1);

// 		Flagt[0 ][j][k].push_back(0);
// 		Flagt[Nx][j][k].push_back(1);
// 	}
// 	Vector3d x0 (0., 0.5*Ny, 0.5*Nz);
// 	Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
// 	Lp[0]->ID = 0;
// 	Lp[0]->Type = -1;
// 	Lp[0]->Normal << 1.,0.,0.;
// 	x0 << Nx, 0.5*Ny, 0.5*Nz;
// 	Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
// 	Lp[1]->ID = 1;
// 	Lp[1]->Type = -1;
// 	Lp[1]->Normal << -1.,0.,0.;

// 	Lp[0]->Lcid.clear();
// 	Lp[1]->Lcid.clear();
// 	for (int j=0; j<=Ny+2; ++j)
// 	for (int k=0; k<=Nz+2; ++k)
// 	{
// 		size_t cid0 = 1+j*Ncy+k*Ncz;
// 		size_t cid1 = Nx+1+j*Ncy+k*Ncz;
// 		Lp[0]->Lcid.push_back(cid0);
// 		Lp[1]->Lcid.push_back(cid1);
// 	}
// 	// For y axis
// 	for (int i=0; i<=Nx; ++i)
// 	for (int k=0; k<=Nz; ++k)
// 	{
// 		Flag[i][0 ][k].push_back(2);
// 		Flag[i][Ny][k].push_back(3);

// 		Flagt[i][0 ][k].push_back(2);
// 		Flagt[i][Ny][k].push_back(3);
// 	}
// 	x0 << 0.5*Nx, 0., 0.5*Nz;		
// 	Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
// 	Lp[2]->ID = 2;
// 	Lp[2]->Type = -1;
// 	Lp[2]->Normal << 0.,1.,0.;
// 	x0 << 0.5*Nx, Ny, 0.5*Nz;
// 	Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
// 	Lp[3]->ID = 3;
// 	Lp[3]->Type = -1;
// 	Lp[3]->Normal << 0.,-1.,0.;

// 	Lp[2]->Lcid.clear();
// 	Lp[3]->Lcid.clear();
// 	for (int i=0; i<=Nx+2; ++i)
// 	for (int k=0; k<=Nz+2; ++k)
// 	{
// 		size_t cid0 = i+Ncy+k*Ncz;
// 		size_t cid1 = i+(Ny+1)*Ncy+k*Ncz;
// 		Lp[2]->Lcid.push_back(cid0);
// 		Lp[3]->Lcid.push_back(cid1);
// 	}

// 	if (D==3)
// 	{
// 		// For z axis
// 		for (int i=0; i<=Nx; ++i)
// 		for (int j=0; j<=Ny; ++j)
// 		{
// 			Flag[i][j][0 ].push_back(4);
// 			Flag[i][j][Nz].push_back(5);

// 			Flagt[i][j][0 ].push_back(4);
// 			Flagt[i][j][Nz].push_back(5);
// 		}
// 	}
// 	x0 << 0.5*Nx, 0.5*Ny, 0.;
// 	Lp.push_back(new DEM_PARTICLE(-5, x0, 0.));
// 	Lp[4]->ID = 4;
// 	Lp[4]->Type = -1;
// 	Lp[4]->Normal << 0.,0.,1.;
// 	x0 << 0.5*Nx, 0.5*Ny, Nz;
// 	Lp.push_back(new DEM_PARTICLE(-6, x0, 0.));
// 	Lp[5]->ID = 5;
// 	Lp[5]->Type = -1;
// 	Lp[5]->Normal << 0.,0.,-1.;

// 	Lp[4]->Lcid.clear();
// 	Lp[5]->Lcid.clear();
// 	for (int i=0; i<=Nx+2; ++i)
// 	for (int j=0; j<=Ny+2; ++j)
// 	{
// 		size_t cid0 = i+j*Ncy+Ncz;
// 		size_t cid1 = i+j*Ncy+(Nz+1)*Ncz;
// 		Lp[4]->Lcid.push_back(cid0);
// 		Lp[5]->Lcid.push_back(cid1);
// 		if (cid0==157)
// 		{
// 			cout <<"add 157 " << endl;
// 			cout << "i= " <<i << " j= " << j << endl;
// 		}
// 	}

// 	MMap.clear();
// 	MMapt.resize(Nproc);
// 	Lsd.resize(Nproc);
// 	cout << "================ Finish init. ================" << endl;
// }

inline void DEM::Init(Vector3d& lx)
{
	Lx = lx;
	int x = (int) ((double) Nx)/Lx(0)+1.;
	int y = (int) ((double) Ny)/Lx(1)+1.;
	int z = (int) ((double) Nz)/Lx(2)+1.;
	Ncz = (x+3)*(y+3);
	Ncy = (x+3);
	Bins.resize((x+3)*(y+3)*(z+3));
	if (D==3)
	{
		NeiV.resize(27);
		for (size_t i=0; i<NeiV.size(); ++i)	NeiV[i] = Nei[i](0)+Nei[i](1)*Ncy+Nei[i](2)*Ncz;	
	}
	Vector3d x0 (0., 0.5*Ny, 0.5*Nz);
	Vector3d n0 (1.,0.,0.);
	AddPlane(-1, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	x0 << Nx, 0.5*Ny, 0.5*Nz;
	n0 << -1., 0., 0.;
	AddPlane(-2, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	x0 << 0.5*Nx, 0., 0.5*Nz;
	n0 << 0., 1., 0.;
	AddPlane(-3, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	x0 << 0.5*Nx, Ny, 0.5*Nz;
	n0 << 0., -1., 0.;
	AddPlane(-4, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	x0 << 0.5*Nx, 0.5*Ny, 0.;
	n0 << 0., 0., 1.;
	AddPlane(-5, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	x0 << 0.5*Nx, 0.5*Ny, Nz;
	n0 << 0., 0., -1.;
	AddPlane(-6, x0, n0);
	Lp[Lp.size()-1]->MID = 1;

	for (int j=0; j<=y+2; ++j)
	for (int k=0; k<=z+2; ++k)
	{
		size_t cid0 = 1+j*Ncy+k*Ncz;
		size_t cid1 = x+j*Ncy+k*Ncz;
		Bins[cid0].push_back(0);
		Bins[cid1].push_back(1);
		// Lp[0]->Lcid.push_back(cid0);
		// Lp[1]->Lcid.push_back(cid1);
	}
	for (int i=0; i<=x+2; ++i)
	for (int k=0; k<=z+2; ++k)
	{
		size_t cid0 = i+Ncy+k*Ncz;
		size_t cid1 = i+y*Ncy+k*Ncz;
		Bins[cid0].push_back(2);
		Bins[cid1].push_back(3);
		// Lp[2]->Lcid.push_back(cid0);
		// Lp[3]->Lcid.push_back(cid1);
	}
	for (size_t i=0; i<=x+2; ++i)
	for (size_t j=0; j<=y+2; ++j)
	{
		size_t cid0 = i+j*Ncy+Ncz;
		size_t cid1 = i+j*Ncy+z*Ncz;
		Bins[cid0].push_back(4);
		Bins[cid1].push_back(5);
		// Lp[4]->Lcid.push_back(cid0);
		// Lp[5]->Lcid.push_back(cid1);
	}
	MMap.clear();
	MMapt.resize(Nproc);
	Lsd.resize(Nproc);
	BinsProc.resize(Nproc);
}

inline void DEM::AddPlane(int tag, Vector3d& x, Vector3d& n)
{
	Lp.push_back(new DEM_PARTICLE(tag, x, 1));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
	Lp[Lp.size()-1]->SetPlane(n);
	Lp[Lp.size()-1]->Fix();
	// UpdateFlag(Lp[Lp.size()-1]);
}

inline void DEM::AddSphere(int tag, double r, Vector3d& x, double rho)
{
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	Lp[p]->SetSphere(r);
	size_t i = ceil(Lp[p]->X(0)/Lx(0))+1;
	size_t j = ceil(Lp[p]->X(1)/Lx(1))+1;
	size_t k = ceil(Lp[p]->X(2)/Lx(2))+1;
	Lp[p]->CID = i+j*Ncy+k*Ncz;		// bin id
	Bins[Lp[p]->CID].push_back(p);
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

inline void DEM::AddMetaball(int tag, double rs, vector<Vector3d>& metaP, VectorXd& metaK, double rho)
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		if (tag==Lp[i]->Tag)
		{
			cout << "\033[1;31mError: Cannot add metaball. The tag is existed.\033[0m\n";		
			exit(0);
		}
	}
	Vector3d x (0,0,0);
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	// cout << Lp.size() << endl;
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	// cout << "SetMetaball" << endl;
	Lp[p]->SetMetaball(rs, metaP, metaK, UnitSpherePoints, UnitSphereFaces);
	size_t i = ceil(Lp[p]->X(0)/Lx(0))+1;
	size_t j = ceil(Lp[p]->X(1)/Lx(1))+1;
	size_t k = ceil(Lp[p]->X(2)/Lx(2))+1;
	Lp[p]->CID = i+j*Ncy+k*Ncz;		// bin id
	Bins[Lp[p]->CID].push_back(p);
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

void DEM::AddNMetaballs(int tag, size_t seed, int np, Vector3d& x0, Vector3d& x1, double rho, double rs, vector<Vector3d>& metaP, VectorXd& metaK)
{
    srand (seed);
    int count = 0;

    Vector3d l = x1-x0;
    cout << l.transpose() << endl;

    for (size_t s=0; s<np; ++s)
    {
    	AddMetaball(Lp.size()+1, rs, metaP, metaK, rho);
    	size_t p = Lp.size()-1;	// ID
    	double r = Lp[p]->R;	// Radius

    	for (int n=0; n<1.0e300; ++n)
    	{
    		bool contacted = false;
	        Vector3d pos = x0;
	        for (size_t d=0; d<D; ++d)	pos(d) += (l(d)-2*(r+1))*((double)rand()/RAND_MAX) + (r+1);
	    	// bin id
	    	size_t cid = ceil(pos(0)/Lx(0))+1+(ceil(pos(1)/Lx(1))+1)*Ncy+(ceil(pos(2)/Lx(2))+1)*Ncz;
			// check overlapping
			for (size_t i=0; i<NeiV.size(); ++i)
			{
				double delta = -100.;
				int b = cid+NeiV[i];			// neighor bins
				for (size_t j=0; j<Bins[b].size(); ++j)
				{
					size_t q = Bins[b][j];
					double delta = -1.;
					if (q<6)	delta = Lp[p]->R-(pos-Lp[q]->X).dot(Lp[q]->Normal);
					else 		delta = Lp[p]->R + Lp[q]->R - (pos-Lp[q]->X).norm();
					if (delta>0.)	contacted = true;
				}
			}
			if (!contacted)
			{
				Vector3d axe;
				axe(0) = (double)rand()/RAND_MAX;
				axe(1) = (double)rand()/RAND_MAX;
				axe(2) = (double)rand()/RAND_MAX;
				axe.normalize();
				double ang = ((double)rand()/RAND_MAX)*2.*M_PI;

				// ang = 0.;
				// pos.setZero();
				// cout << "ang: " << ang << endl;
				// cout << "axe: " << axe.transpose() << endl;
				Matrix3d m;
				m = AngleAxisd(ang, axe);
				Quaterniond q0(m);
				Lp[p]->Q0 = q0;
				Lp[p]->X = pos;
				for (size_t i=0; i<Lp[p]->MetaP.size(); ++i)
				{
					Lp[p]->MetaP[i] = q0._transformVector(Lp[p]->MetaP0[i]);
					Lp[p]->MetaP[i] += pos;
				}
				for (size_t i=0; i<Lp[p]->P.size(); ++i)
				{
					Lp[p]->P[i] = q0._transformVector(Lp[p]->P0[i]);
					Lp[p]->P[i] += pos;
				}
				Lp[p]->CID = cid;		// bin id
				Bins[cid].push_back(p);
				// cout << "cid: " << cid << endl;
				cout << "add metaball: " << Lp[p]->Tag << endl;
				break;
			}
    	}
    }
}

inline void DEM::Move(bool movePoints)
{
	// #pragma omp parallel for schedule(dynamic, 1000) num_threads(Nproc)
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
			p0->VelocityVerlet(Dt, movePoints);
			p0->UpdateBox(D);
			// if (p0->V.norm()>0.003)
			// if (p0->X(2)<0.5)
			// {
			// 	cout << "p0->ID: " << p0->ID << endl;
			// 	cout << "p0->CID: " << p0->CID << endl;
			// 	abort();
			// }
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
		g0->VelocityVerlet(Dt, movePoints);
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

inline void DEM::UpdatePositionGlobal(double dt, bool movePoints, bool resetHydro)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		DEM_PARTICLE* p0 = Lp[i];
		p0->UpdatePosition(dt, movePoints, resetHydro);
	}
}

inline void DEM::UpdateVelocityGlobal(double dt)
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=6; i<Lp.size(); ++i)
	{
		DEM_PARTICLE* p0 = Lp[i];
		p0->UpdateVelocity(dt);
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
	// gn 	= -2.*Beta*sqrt(kn*me);
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
inline void DEM::Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, Vector3d& xir, Vector3d& finalP, bool& contacted)
{
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
	if ((pi->Type==4 && pj->Type==4) || (pi->Type==1 && pj->Type==4) || (pi->Type==4 && pj->Type==1))			// For collision of metaballs
	{
		if (pi->Type==1)
		{
			pi->MetaP[0] = pi->X;
		}
		else if (pj->Type==1)
		{
			pj->MetaP[0] = pj->X;
		}
		bool conv = true;
		bool useCell = false;
		Vector3d initP;						// inital point for newton's method
		size_t key = Key(pi->ID, pj->ID);
		if (MMap.count(key)>0)
		{
			initP = MMap[key];
			FindMetaClosestPoints(pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, initP, finalP, cpi, cpj, conv);
		}// use existing one
		else 										// for new pair find closest control point pair
		{
			// size_t ii,jj;
			// double minDis2 = 1.0e30;
			// for (size_t i=0; i<pi->MetaP.size(); ++i)
			// for (size_t j=i; j<pj->MetaP.size(); ++j)
			// {
			// 	double dis2 = (pi->MetaP[i]-pj->MetaP[j]).squaredNorm();
			// 	if (dis2<minDis2)
			// 	{
			// 		ii = i;
			// 		jj = j;
			// 		minDis2 = dis2;
			// 	}
			// }
			FindMetaClosestPointsWithCell(Ne, pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, 8, finalP, cpi, cpj);
			useCell = true;
			// initP = 0.5*(pi->MetaP[ii]+pj->MetaP[jj]);
		}
		// FindMetaClosestPoints(pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, initP, finalP, cpi, cpj, conv);
		if (!conv)		// if not converge
		{
			cout << "pi->ID: " << pi->ID << endl;
			cout << "pj->ID: " << pj->ID << endl;
			FindMetaClosestPointsWithCell(Ne, pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, 8, finalP, cpi, cpj);
			useCell = true;
			abort();
			// use cell method to find inital point
			// FindMetaInitPoint(Ne, pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, 5, initP);
			// conv = true;
			// FindMetaClosestPoints(pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, initP, finalP, cpi, cpj, conv);
			// if (!conv)
			// {
			// 	// FindMetaInitPoint(Ne, pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, 8, initP);
			// 	// cout << "still not converge!" << endl;
			// 	// for (size_t i=0; i<pi->MetaP.size(); ++i)
			// 	// {
			// 	// 	cout << pi->MetaP[i].transpose() << endl;
			// 	// 	cout << "ki: " << pi->MetaK[i] << endl;
			// 	// 	cout << pj->MetaP[i].transpose() << endl;
			// 	// 	cout << "kj: " << pj->MetaK[i] << endl;
			// 	// }
			// 	FindMetaClosestPointsWithCell(Ne, pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, 8, finalP, cpi, cpj);
			// 	useCell = true;
			// 	// abort();
			// 	// conv = true;
			// 	// FindMetaClosestPoints(pi->MetaP, pj->MetaP, pi->MetaK, pj->MetaK, initP, MMapt[tid][key], cpi, cpj, conv);
			// }
		}
		// int tid = omp_get_thread_num();
		// MMapt[tid][key] = finalP;
		n = cpi-cpj;
		delta = pi->Rs+pj->Rs-n.norm();
		n.normalize();
		cp = 0.5*(cpi+cpj);
		// if (delta>0.)	contacted = true;
		// cout << "delta: " << delta << endl;
		// if (delta>-0.0001)
		// if (pi->Tag ==120 && pj->Tag==284 && delta>0.)
		// if (pi->ID ==151 && pj->ID==155 /*&& delta>0.*/)
		// {
		// 	if (useCell)	cout << "useCell" << endl;
		// 	// cout << "pi->ID: " << pi->ID << endl;
		// 	// cout << "pj->CID: " << pj->CID << endl;
		// 	cout << pi->R+pj->R-(pi->X-pj->X).norm() << endl;
		// 	cout << "delta: " << delta << endl;
			// cout << "cpi: " << cpi.transpose() << endl;
		// 	// double ci = CalMetaC(pi->MetaP, pi->MetaK, cpi);
		// 	double ci, cj;
		// 	Vector3d fi, fj;
		// 	CalMetaCF(pi->MetaP, pi->MetaK, cpi, ci, fi);
		// 	CalMetaCF(pj->MetaP, pj->MetaK, cpj, cj, fj);
			// cout << "ci: " << CalMetaC(pi->MetaP, pi->MetaK, cpi) << endl;
			// cout << "cpj: " << cpj.transpose() << endl;
			// cout << "cj: " << CalMetaC(pj->MetaP, pj->MetaK, cpj) << endl;
		// 	cout << "n: " << n.transpose() << endl;
		// 	for (size_t s=0; s<pi->MetaP.size(); ++s)
		// 	{
		// 		cout << "pi->MetaP: " << pi->MetaP[s].transpose() << endl;
		// 		cout << "pj->MetaP: " << pj->MetaP[s].transpose() << endl;
		// 		cout << "pi->MetaK: " << pi->MetaK[s] << endl;
		// 	}
		// 	cout << "initP: " << initP.transpose() << endl;
		// 	cout << "finalP: " << finalP.transpose() << endl;
		// 	cout << "pi->V: " << pi->V.transpose() << endl;
		// 	cout << "pj->V: " << pj->V.transpose() << endl;
		// 	// if (ci<0.5)	abort();
		// 	cout << "cpi= " << cpi.transpose() << endl;
		// 	cout << "cpj= " << cpj.transpose() << endl;
		// 	n = cpi-cpj;
		// 	n.normalize();
		// 	cout << "n= " << n.transpose() << endl;
		// 	abort();
		// }
	}
	else if (pi->Type==-1 && pj->Type==4)			// For collision of metaball and plane
	{
		size_t key = Key(pi->ID, pj->ID);
		bool first = false;
		if (MMap.count(key)==0)	first = true;
		bool conv = true;
		FindMetaPlaneClosestPoints(pi->Normal, pi->X, pj->MetaP, pj->MetaK, pj->X, first, MMap[key], finalP, cpi, cpj, conv);
		if (!conv)
		{
			FindMetaPlaneClosestPointsWithCell(pi->Normal, pi->X, pj->MetaP, pj->MetaK, pj->X, 5, finalP, cpi, cpj);
			// cout << "using FindMetaPlaneClosestPointsWithCell" << endl; 
			// first = false;
			// FindMetaPlaneClosestPoints(pi->Normal, pi->X, pj->MetaP, pj->MetaK, pj->X, first, initP, MMapt[tid][key], cpi, cpj, conv);
			// abort();
		}
		// int tid = omp_get_thread_num();
		// MMapt[tid][key] = finalP;
		n = cpi-cpj;
		delta = pj->Rs-n.norm();

		n.normalize();
		cp = 0.5*(cpi+cpj);
		// if (delta>0.)	contacted = true;

		if (!conv)
		{
			cout << "pi->Tag: " << pi->Tag << endl;
			cout << "pj->Tag: " << pj->Tag << endl;
			cout << "pi->ID: " << pi->ID << endl;
			cout << "pj->ID: " << pj->ID << endl;
			cout << "delta: " << delta << endl;	
			cout << "cpi: " << cpi.transpose() << endl;
			cout << "cpj: " << cpj.transpose() << endl;		
			for (size_t s=0; s<pj->MetaP.size(); ++s)
			{
				// cout << "pi->MetaP: " << pi->MetaP[s].transpose() << endl;
				cout << "pj->MetaP: " << pj->MetaP[s].transpose() << endl;
				cout << "pj->MetaK: " << pj->MetaK[s] << endl;
			}
			abort();
		}
		// if (delta>0. && pj->Tag==4)
		// {
		// 	cout << "detla: " << delta << endl;
		// }
		// if (delta>-0.001 && pj->Tag==6)
		// {
		// 	cout << "detla 4: " << delta << endl;
		// 	cout << "n: " << n.transpose() << endl;
		// }
		// if (pj->Tag==1)
		// {
		// 	cout << "detla: " << delta << endl;
		// 	if (delta<tmpd && !conv /*&& delta<-0.32*/)	abort();
		// 	tmpd = delta;
		// 	// abort();
		// }

		// if (pi->ID ==31 && pj->ID==63 /*&& delta>0.*/)
		// {
		// 	// cout << "pi->ID: " << pi->ID << endl;
		// 	// cout << "pj->CID: " << pj->CID << endl;
		// 	cout << "delta: " << delta << endl;
		// 	// abort();
		// }
		// if (pi->ID ==3 && pj->ID==85 /*&& delta>0.*/)
		// {
		// 	// if (useCell)	cout << "useCell" << endl;
		// 	// cout << "pi->ID: " << pi->ID << endl;
		// 	// cout << "pj->CID: " << pj->CID << endl;
		// 	cout << pi->R+pj->R-(pi->X-pj->X).norm() << endl;
		// 	cout << "delta: " << delta << endl;
		// 	cout << "cpi: " << cpi.transpose() << endl;
		// 	cout << "cpj: " << cpj.transpose() << endl;
		// 	cout << "n: " << n.transpose() << endl;
		// 	// for (size_t s=0; s<pi->MetaP.size(); ++s)
		// 	// {
		// 	// 	// cout << "pi->MetaP: " << pi->MetaP[s].transpose() << endl;
		// 	// 	cout << "pj->MetaP: " << pj->MetaP[s].transpose() << endl;
		// 	// 	// cout << "pi->MetaK: " << pi->MetaK[s] << endl;
		// 	// }
		// 	// cout << "initP: " << initP.transpose() << endl;
		// 	cout << "finalP: " << finalP.transpose() << endl;
		// 	// abort();
		// }
	}
	// else if (pi->Type==4 && pj->Type==-1)
	// {
	// 	FindMetaPlanceClosestPoints(pj->Normal, pj->X, pi->MetaP, pi->MetaK, pi->X, cpj, cpi);
	// }

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

		// if (delta>0.02)	cout << "delta larger than Rs: " << delta << endl;
// if (pj->Tag==4)
// {		cout << "delta: " << delta << endl;
// 		cout << "kn: " << kn << endl;
// 		cout << "gn: " << gn << endl;
// 		cout << "fn: " << fn.transpose() << endl;
// 		cout << "n: " << n.transpose() << endl;
// 		cout << "vn: " << vn.transpose() << endl;
// 		cout << "pi->X: " << pi->X.transpose() << endl;
// 		cout << "pj->X: " << pj->X.transpose() << endl;
// 		cout << "pi->V: " << pi->V.transpose() << endl;
// 		cout << "pj->V: " << pj->V.transpose() << endl;}

// 		cout << "fn/m: " << fn.transpose()/pi->M << endl;
// 		cout << "m: " << pi->M << endl;
// 		cout << "Vol: " << pi->Vol << endl;
// 		cout << "Rho: " << pi->Rho << endl;
		// abort();


		Vector3d ft (0., 0., 0.);
		// Friction(pi, pj, delta, kt, gt, n, fn, xi, ft);	// Friction force and torque
		Friction(pi, pj, kt, gt, cp, n, fn, xi, ft);
		// Rolling resistance torque
		// double kr = RatioKnr*kn;
		// double gr = RatioGnr*gn;
		// RollingResistance(pi, pj, delta, kr, gr, n, fn, xir, armr);
		// arm += armr;
		Vector3d fnt = fn+ft;							// Total force
		Vector3d tci = (pi->Qfi._transformVector(fnt)).cross(pi->Qfi._transformVector(Xi-cp));
		Vector3d tcj = (pj->Qfi._transformVector(fnt)).cross(pj->Qfi._transformVector(Xj-cp));

		// fnt = fn;
		// tci.setZero();
		// tcj.setZero();
		
		for (size_t d=0; d<D; ++d)
		{
			#pragma omp atomic
			pi->Fc(d) += fnt(d);
			#pragma omp atomic
			pj->Fc(d) -= fnt(d);
			#pragma omp atomic
			pi->Tc(d) += tci(d);
			#pragma omp atomic
			pj->Tc(d) -= tcj(d);
		}	
	}
}

inline void DEM::Friction(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double kt, double gt, Vector3d& cp, Vector3d& n, Vector3d& fn, Vector3d& xi, Vector3d& ft)
{
	// Relative velocity at the contact point
	Vector3d wi = pi->Qf._transformVector(pi->W);
	Vector3d wj = pj->Qf._transformVector(pj->W);

	Vector3d vij = pi->V + wi.cross(cp-pi->X) - pj->V - wj.cross(cp-pj->X);

	// Vector3d vij = pi->V-pj->V+(pi->R-0.5*delta)*n.cross(pi->W)+(pj->R-0.5*delta)*n.cross(pj->W);	// eq.10
	// Vector3d vij = pi->V-pj->V+(pi->R-0.5*delta)*n.cross(pi->W)+(pj->R-0.5*delta)*n.cross(pj->W);	// eq.10
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
		// else if (p0->Type==4)	dis =p0->MetaCC-CalMetaC(p0->MetaP, p0->MetaK, ind);
		else if (p0->Type==4)	dis = (ind-p0->X).norm()-p0->R-0.87;
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

// inline void DEM::LinkedCell()
// {
// 	CMap.clear();
// 	// sort particles into processor bins for omp
// 	vector<vector<size_t>> lps;
// 	for (size_t p=6; p<Lp.size(); ++p)
// 	{
// 		lps[Lp[p]->PID].push_back(p);
// 	}
// 	vector<vector <vector<size_t>>> lcProc(Nproc);	// contact list for each processor
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t c=0; c<Nproc; ++c)
// 	{
// 		// copy info to local, maybe faster
// 		size_t ncz = Lsd[c]->Ncz;
// 		size_t ncy = Lsd[c]->Ncy;
// 		int nc = Lsd[c]->Nc;
// 		Vector3d xmin = Lsd[c]->Xmin;
// 		Vector3d xmax = Lsd[c]->Xmax;
// 		Vector3d lx = Lx;

// 		vector<vector<size_t>> cells; 	// List of cells
// 		unordered_map<size_t, bool> cmap;
// 		cells.resize(nc);
// 		// sort particles to cells
// 		for (size_t n=0; n<lps[c].size(); ++n)
// 		{
// 			size_t p = lps[c][n];
// 			Vector3d xr = Lp[p]->X - xmin;	// relative position 
// 			Vector3i ind;					// index of the 
// 			ind(0) = (int) xr(0)/lx(0);
// 			ind(1) = (int) xr(1)/lx(1);
// 			ind(2) = (int) xr(2)/lx(2);
// 			int cid = ind(0)+ind(1)*ncy+ind(2)*ncz;	// 
// 			cells[cid].push_back(p);

// 			Vector3i indn = ind;
// 			vector<bool> ord(0);
// 			int sumOrd = 0;
// 			for (size_t d=0; d<D; ++d)
// 			{
// 				bool border = true;
// 				if (xr(d)<Lp[p]->R)					indn(d) -= 1;
// 				else if (xr(d)+Lp[p]->R>xmax(d))	indn(d) += 1;
// 				else 								border   = 0;
// 				ord.push_back(border);
// 				sumOrd += border;
// 			}
// 			if (sumOrd==1)
// 			{
// 				size_t cid1 = indn(0)+indn(1)*ncy+indn(2)*ncz;	// face
// 				cells[cid1].push_back(p);
// 			}
// 			else if (sumOrd==2)
// 			{
// 				int cid1 = ind(0)+indn(1)*ncy+indn(2)*ncz;	// edge
// 				int cid2 = indn(0)+ind(1)*ncy+indn(2)*ncz;	// edge
// 				int cid3 = indn(0)+indn(1)*ncy+ind(2)*ncz;	// edge
// 				if (cid1>=0 && cid1<=nc) 	cells[cid1].push_back(p);
// 				if (cid2>=0 && cid2<=nc)	cells[cid2].push_back(p);
// 				if (cid3>=0 && cid3<=nc)	cells[cid3].push_back(p);
// 			}
// 			else
// 			{
// 				int cid0 = indn(0)+indn(1)*ncy+indn(2)*ncz;	// corner
// 				int cid1 = ind(0)+indn(1)*ncy+indn(2)*ncz;	// edge
// 				int cid2 = indn(0)+ind(1)*ncy+indn(2)*ncz;	// edge
// 				int cid3 = indn(0)+indn(1)*ncy+ind(2)*ncz;	// edge	
// 				int cid4 = ind(0)+ind(1)*ncy+indn(2)*ncz;	// face
// 				int cid5 = indn(0)+ind(1)*ncy+ind(2)*ncz;	// face
// 				int cid6 = ind(0)+indn(1)*ncy+ind(2)*ncz;	// face

// 				if (cid0>=0 && cid0<=nc) 	cells[cid0].push_back(p);
// 				if (cid1>=0 && cid1<=nc)	cells[cid1].push_back(p);
// 				if (cid2>=0 && cid2<=nc)	cells[cid2].push_back(p);
// 				if (cid3>=0 && cid3<=nc)	cells[cid3].push_back(p);
// 				if (cid4>=0 && cid4<=nc)	cells[cid4].push_back(p);
// 				if (cid5>=0 && cid5<=nc)	cells[cid5].push_back(p);
// 				if (cid6>=0 && cid6<=nc)	cells[cid6].push_back(p);
// 			}
// 		}
// 		// find collision from every cell
// 		for (size_t s=0; s<cells.size(); ++s)
// 		{
// 			for (size_t i=0; i<cells[s].size(); ++i)
// 			for (size_t j=i+1; j<cells[s].size(); ++j)
// 			{
// 				size_t p = min(cells[s][i],cells[s][j]);
// 				size_t q = max(cells[s][i],cells[s][j]);
// 				double delta;
// 				if (Lp[p]->Type==-1 && Lp[q]->Type!=-1)		delta = Lp[q]->R - (Lp[q]->X-Lp[p]->X).dot(Lp[p]->Normal);
// 				else										delta = Lp[p]->R+Lp[q]->R - (Lp[p]->X-Lp[q]->X).norm();
// 				size_t key = Key(p, q);
// 				if (delta>0. && (!cmap[key]))
// 				{
// 					lcProc[c].push_back({p, q});
// 					cmap[key] = true;
// 				}
// 			}
// 		}
// 	}
// 	// combine contact lists together
// 	for (size_t c=0; c<Nproc; ++c)
// 	for (size_t i=0; i<lcProc[c].size(); ++i)
// 	{
// 		// size_t p = lcProc[c][0];
// 		// size_t p = lcProc[c][1];
// 		size_t key = Key(lcProc[c][i][0], lcProc[c][i][1]);
// 		if (!CMap[key])
// 		{
// 			Lc.push_back(lcProc[c][i]);
// 			CMap[key] = true;
// 		}
// 	}
// }

// inline void DEM::LinkedCell()
// {
// 	CMap.clear();
// 	// 3d data structure of cells, first index: proc, second: bins, third: particle id
// 	// the idea is to make sure every processor have their own bins
// 	vector<vector<vector<size_t>>> binProc(Nproc);
// 	vector<unordered_map<size_t,size_t>>	bmapProc(Nproc);
// 	for (size_t c=0; c<Nproc; ++c)	bmapProc[c].clear();

// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=6; p<Lp.size(); ++p)
// 	{
// 		size_t tid = omp_get_thread_num();				// thread id
// 		Vector3d xp = Lp[p]->X;
// 		// Vector3d xr (xp(0)/Lx(0), xp(1)/Lx(1), xp(2)/Lx(2));
// 		// Vector3d rr (Lp[p]->R/Lx(0), Lp[p]->R/Lx(1), Lp[p]->R/Lx(2));
// 		Vector3i id;									// bin id
// 		id(0) = (int) xr(0);
// 		id(1) = (int) xr(1);
// 		id(2) = (int) xr(2);

// 		// vector<int> lcid(0);							// list of bin id that need to push p into it
// 		int cid = id(0)+id(1)*Ncy+id(2)*Ncz;
// 		// lcid.push_back(cid);
// 		if (bmapProc[tid].count(cid)==0)				// if this bin is new for this thread
// 		{
// 			bmapProc[tid][cid] = binProc[tid].size();	// mark this bin in bin map
// 			binProc[tid].push_back({p});				// add to this thread's bin system
// 		}
// 		else 											// if this bin already exist
// 		{
// 			size_t bid = bmapProc[tid][cid];			// find the id in bin system
// 			binProc[tid][bid].push_back(p);				// add this particle to this bin
// 		}	
// 		// Vector3i idn = id;
// 		// vector<bool> ord(0);
// 		// int sumOrd = 0;
// 		// for (size_t d=0; d<D; ++d)
// 		// {
// 		// 	bool border = true;
// 		// 	if (xr(d)-id(d)<rr(d))			idn(d) -= 1;
// 		// 	else if (xr(d)-id(d)>1.-rr(d))	idn(d) += 1;
// 		// 	else 							border  = 0;
// 		// 	ord.push_back(border);
// 		// 	sumOrd += border;
// 		// }

// 		// if (sumOrd==1)
// 		// {
// 		// 	size_t cid1 = idn(0)+idn(1)*Ncy+idn(2)*Ncz;	// face
// 		// 	lcid.push_back(cid1);
// 		// }
// 		// else if (sumOrd==2)
// 		// {
// 		// 	int cid1 = id(0)+idn(1)*Ncy+idn(2)*Ncz;	// edge
// 		// 	int cid2 = idn(0)+id(1)*Ncy+idn(2)*Ncz;	// edge
// 		// 	int cid3 = idn(0)+idn(1)*Ncy+id(2)*Ncz;	// edge
// 		// 	if (cid1>=0 && cid1<=Nc) 	lcid.push_back(cid1);
// 		// 	if (cid2>=0 && cid2<=Nc)	lcid.push_back(cid2);
// 		// 	if (cid3>=0 && cid3<=Nc)	lcid.push_back(cid3);
// 		// }
// 		// else
// 		// {
// 		// 	int cid0 = idn(0)+idn(1)*Ncy+idn(2)*Ncz;// corner
// 		// 	int cid1 = id(0)+idn(1)*Ncy+idn(2)*Ncz;	// edge
// 		// 	int cid2 = idn(0)+id(1)*Ncy+idn(2)*Ncz;	// edge
// 		// 	int cid3 = idn(0)+idn(1)*Ncy+id(2)*Ncz;	// edge	
// 		// 	int cid4 = id(0)+id(1)*Ncy+idn(2)*Ncz;	// face
// 		// 	int cid5 = idn(0)+id(1)*Ncy+id(2)*Ncz;	// face
// 		// 	int cid6 = id(0)+idn(1)*Ncy+id(2)*Ncz;	// face

// 		// 	if (cid0>=0 && cid0<=Nc) 	lcid.push_back(cid0);
// 		// 	if (cid1>=0 && cid1<=Nc)	lcid.push_back(cid1);
// 		// 	if (cid2>=0 && cid2<=Nc)	lcid.push_back(cid2);
// 		// 	if (cid3>=0 && cid3<=Nc)	lcid.push_back(cid3);
// 		// 	if (cid4>=0 && cid4<=Nc)	lcid.push_back(cid4);
// 		// 	if (cid5>=0 && cid5<=Nc)	lcid.push_back(cid5);
// 		// 	if (cid6>=0 && cid6<=Nc)	lcid.push_back(cid6);
// 		// }

// 		// for (size_t i=0;i<lcid.size(); ++i)
// 		// {
// 		// 	if (bmapProc[tid].count(lcid[i])==0)				// if this bin is new for this thread
// 		// 	{
// 		// 		bmapProc[tid][lcid[i]] = binProc[tid].size();	// mark this bin in bin map
// 		// 		binProc[tid].push_back({p});					// add to this thread's bin system
// 		// 	}
// 		// 	else 												// if this bin already exist
// 		// 	{
// 		// 		size_t bid = bmapProc[tid][lcid[i]];			// find the id in bin system
// 		// 		binProc[tid][bid].push_back(p);					// add this particle to this bin
// 		// 	}			
// 		// }
// 	}

// 	#pragma omp parallel for schedule(dynamic,4) num_threads(Nproc)
// 	for (size_t i=0; i<Bins.size(); ++i)
// 	{
// 		Bins[i].clear();
// 		for (size_t c=0; c<Nproc; ++c)
// 		{
// 			if (bmapProc[c].count(i)!=0)
// 			{
// 				size_t bid = bmapProc[c][i];
// 				Bins[i].insert(binProc[c][bid].begin(),binProc[c][bid].end());
// 			}
// 		}
// 	}

// 	#pragma omp parallel for schedule(dynamic,4) num_threads(Nproc)
// 	for (size_t i=0; i<Bins.size(); ++i)
// 	{

// 	}
// }

// inline void DEM::SetLinkedCell(Vector3d& lx)
// {
// 	Lx = lx;
// 	int x = (int) ((double) Nx)/Lx(0)+1.;
// 	int y = (int) ((double) Ny)/Lx(1)+1.;
// 	int z = (int) ((double) Nz)/Lx(2)+1.;
// 	Ncz = (x+3)*(y+3);
// 	Ncy = (x+3);
// 	Bins.resize((x+3)*(y+3)*(z+3));
// 	if (D==3)
// 	{
// 		NeiV.resize(27);
// 		for (size_t i=0; i<NeiV.size(); ++i)	NeiV[i] = Nei[i](0)+Nei[i](1)*Ncy+Nei[i](2)*Ncz;	
// 	}
// 	cout << "nx: " << nx << endl;
// 	cout << "ny: " << ny << endl;
// 	cout << "nz: " << nz << endl;
// 	abort();
// }

// inline void DEM::LinkedCell()
// {
// 	CMap.clear();
// 	// 3d data structure of cells, first index: proc, second: bins, third: particle id
// 	// the idea is to make sure every processor have their own bins
// 	vector<vector<vector<size_t>>> binProc(Nproc);
// 	vector<unordered_map<size_t,size_t>>	bmapProc(Nproc);
// 	for (size_t c=0; c<Nproc; ++c)	bmapProc[c].clear();

// 	vector<vector<vector<size_t>>> 			binsProc(Nproc);
// 	for (size_t c=0; c<Nproc; ++c)	binsProc[c].resize(Bins.size());
// 	// cout << "sizeof(binsProc): " << sizeof(binsProc) << endl;
// 	// cout << "sizeof(Bins): " << sizeof(Bins) << endl;
// 	// abort();

// 	// size_t DomSplit[nx+3][ny+3][nz+3];

// 	// cout << "111111111" << endl;
// 	auto t_start = std::chrono::system_clock::now();
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=6; p<Lp.size(); ++p)
// 	{
// 		size_t tid = omp_get_thread_num();				// thread id
// 		Vector3d xp = Lp[p]->X+Lx;						// move to extended domain to avoid boundary issue
// 		Vector3i ind;									// bin 3d index
// 		ind(0) = (int) xp(0)/Lx(0);
// 		ind(1) = (int) xp(1)/Lx(1);
// 		ind(2) = (int) xp(2)/Lx(2);

// 		int cid = ind(0)+ind(1)*Ncy+ind(2)*Ncz;			// bin id
// 		Lp[p]->CID = cid;
// 		Lp[p]->Lbins.clear();
// 		for (size_t i=0; i<Nei.size(); ++i)
// 		{
// 			Vector3i indi = ind + Nei[i];
// 			int cidi = indi(0)+indi(1)*Ncy+indi(2)*Ncz;
// 			Lp[p]->Lbins.push_back(cidi);
// 		}
// 		// binsProc[tid][cid].push_back(p);
// 		// #pragma omp critical
// 		// {
// 		// 	Bins[cid].push_back(p);
// 		// }
// 		// Bins[cid].push_back(p);
// 		if (bmapProc[tid].count(cid)==0)				// if this bin is new for this thread
// 		{
// 			bmapProc[tid][cid] = binProc[tid].size();	// mark this bin in bin map
// 			binProc[tid].push_back({p});				// add to this thread's bin system
// 		}
// 		else 											// if this bin already exist
// 		{
// 			size_t bid = bmapProc[tid][cid];			// find the id in bin system
// 			binProc[tid][bid].push_back(p);				// add this particle to this bin
// 		}
// 	}
// 	auto t_end = std::chrono::system_clock::now();
// 	// cout << "111 time: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 	// cout << "222222222" << endl;
// 	// combine bins for each processor togther
// 	t_start = std::chrono::system_clock::now();
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t i=0; i<Bins.size(); ++i)
// 	{
// 		Bins[i].clear();
// 		for (size_t c=0; c<Nproc; ++c)
// 		{
// 			if (bmapProc[c].count(i)!=0)				// if this processor have this bin id
// 			{
// 				size_t bid = bmapProc[c][i];			// get bin position
// 				Bins[i].insert(Bins[i].end(),binProc[c][bid].begin(),binProc[c][bid].end());
// 			}
// 		}
// 	}
// 	t_end = std::chrono::system_clock::now();
// 	// cout << "222 time: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 	// // cout << "3333333333" << endl;
// 	t_start = std::chrono::system_clock::now();
// 	vector<vector<vector<size_t>>> lcProc(Nproc);
// 	// loop over particles to check collision
// 	#pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (size_t p=6; p<Lp.size(); ++p)
// 	{
// 		size_t tid = omp_get_thread_num();
// 		for (size_t i=0; i<Lp[p]->Lbins.size(); ++i)
// 		{
// 			// cout << "p: " << p << endl;
// 			int b = Lp[p]->Lbins[i];
// 			// cout << "b: " << b << endl;
// 			for (size_t j=0; j<Bins[b].size(); ++j)
// 			{
// 				// cout << "j: " << j << endl;
// 				// cout << Bins[b].size() << endl;
// 				size_t q = Bins[b][j];
// 				// cout << "q: " << q << endl;
// 				if (q>p)
// 				{
// 					double delta = Lp[p]->R + Lp[q]->R - (Lp[p]->X-Lp[q]->X).norm();
// 					if (delta>0.)	lcProc[tid].push_back({p,q});
// 				}
// 			}
// 		}
// 	}
// 	t_end = std::chrono::system_clock::now();
// 	// cout << "333 time: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 	// cout << "44444444444" << endl;
// 	t_start = std::chrono::system_clock::now();
// 	Lc.clear();
// 	for (size_t c=0; c<Nproc; ++c)	Lc.insert(Lc.end(), lcProc[c].begin(),lcProc[c].end());
// 	t_end = std::chrono::system_clock::now();
// 	// cout << "444 time: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 	// // abort();
// }
inline void DEM::FindIndex(size_t n, size_t& i, size_t& j, size_t& k)
{
	k = n/Ncz;
	j = (n%Ncz)/Ncy;
	i = (n%Ncz)%Ncy;
}

/*inline void DEM::LinkedCell()
{
	CMap.clear();

	vector< vector<size_t> > updateLpProc(Nproc);		// list of particles that need to update it's bin id
	// cout << "111111111" << endl;
	auto t_start = std::chrono::system_clock::now();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		size_t tid = omp_get_thread_num();				// thread id

		size_t i = ceil(Lp[p]->X(0)/Lx(0))+1;
		size_t j = ceil(Lp[p]->X(1)/Lx(1))+1;
		size_t k = ceil(Lp[p]->X(2)/Lx(2))+1;		

		Lp[p]->CIDt = i+j*Ncy+k*Ncz;			// bin id
		if (Lp[p]->CID != Lp[p]->CIDt)	updateLpProc[tid].push_back(p);
	}
	auto t_end = std::chrono::system_clock::now();
	// size_t nm = 0;
	for (size_t c=0; c<Nproc; ++c)
	{
		if (updateLpProc[c].size()>0)
		{
			for (size_t i=0; i<updateLpProc[c].size(); ++i)
			{
				size_t p = updateLpProc[c][i];
				for (size_t j=0; j<Bins[Lp[p]->CID].size(); ++j)
				{
					if (Bins[Lp[p]->CID][j]==p)
					{
						Bins[Lp[p]->CID].erase(Bins[Lp[p]->CID].begin() + j);
						break;
					}
				}
				// cout << Lp[p]->X.transpose() << endl;
				Bins[Lp[p]->CIDt].push_back(p);
				Lp[p]->CID = Lp[p]->CIDt;
				// nm++;
			}			
		}
	}
	vector<vector<vector<size_t>>> lcProc(Nproc);
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		size_t tid = omp_get_thread_num();
		bool walls[6] {false, false, false, false, false, false};
		for (size_t i=0; i<NeiV.size(); ++i)
		{
			int b = Lp[p]->CID+NeiV[i];			// neighor bins
			// if (b==301)	cout << "bin: 301" << endl;
			for (size_t j=0; j<Bins[b].size(); ++j)
			{
				size_t q = Bins[b][j];
				if (q<6)
				{
						// if (q==1 && p==121)
						// {
						// 	cout << "q000: " << q << endl;
						// 	// cout << "delta++++: " << delta << endl;							
						// }
					// cout << "q: " << q << endl;
					if (!walls[q])
					{
						double delta = Lp[p]->R-(Lp[p]->X-Lp[q]->X).dot(Lp[q]->Normal);
						if (delta>0.)
						{
							lcProc[tid].push_back({q,p});
							walls[q] = true;
						}
						// if (q==4 && p==68)
						// {
						// 	cout << "q: " << q << endl;
						// 	cout << "delta++++: " << delta << endl;							
						// }
						// cout << "q: " << q << endl;
						// cout << "delta: " << delta << endl;
					}
				}
				else if (q>p)
				{
					double delta = Lp[p]->R + Lp[q]->R - (Lp[p]->X-Lp[q]->X).norm();
					if (delta>0.)	lcProc[tid].push_back({p,q});
				}
			}
		}
	}
	// cout << "44444" << endl;
	Lc.clear();
	for (size_t c=0; c<Nproc; ++c)	Lc.insert(Lc.end(), lcProc[c].begin(),lcProc[c].end());
		// cout << "555555555" << endl;
	// for (size_t i=0; i<Lp[4]->Lcid.size(); ++i)
	// {
	// 	size_t b = Lp[4]->Lcid[i];
	// 	if (b==157)
	// 	{
	// 		cout << "b=157" << endl;
	// 		// abort();
	// 	}
	// }
}*/

// Generate BinsProc for openMP with balanced load
inline void DEM::GenBinsProc()
{
	size_t ab = 0;		// actived bin number
	size_t abe = 0;		// actived bin element number
	for (size_t b=0; b<Bins.size(); ++b)
	{
		if (Bins[b].size()>0)
		{
			ab++;
			abe += Bins[b].size();
		}
	}
	size_t av = ceil(((double) abe)/((double) Nproc))+1;	// average element number in a proc
	size_t c = 0;
	size_t d = 0;
	for (size_t b=0; b<Bins.size(); ++b)
	{
		if (d>av && c<Nproc)	// stop add if element number larger than average
		{
			c++;
			d = 0;
		}
		if (Bins[b].size()>0)
		{
			BinsProc[c].push_back(b);
			d += Bins[b].size();
		}
	}
	for (size_t c=0; c<Nproc; ++c)
	{
		size_t num = 0;
		for (size_t s=0; s<BinsProc[c].size(); ++s)
		{
			size_t b = BinsProc[c][s];
			num += Bins[b].size();
		}
		// cout << "c: " << c << endl;
		// cout << "BinsProc[c].size(): " << BinsProc[c].size() << endl;
		// cout << "num: " << num << endl;
		// cout << "++++++++++++" << endl;
	}
}

inline void DEM::LinkedCell()
{
	CMap.clear();

	vector< vector<size_t> > updateLpProc(Nproc);		// list of particles that need to update it's bin id
	// Loop for all bins and find particles that need to be update their bin ID
	size_t BinNump[Nproc];
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<Nproc; ++c)
	{
		BinNump[c] = 0;
		for (size_t s=0; s<BinsProc[c].size(); ++s)
		{
			size_t b = BinsProc[c][BinsProc[c].size()-1-s];	// bin id
			if (Bins[b].size()==0)							// remove if empty bin
			{
				swap(BinsProc[c][s],BinsProc[c].back());
				BinsProc[c].pop_back();
			}
			else
			{
				BinNump[c] += Bins[b].size();
				// Loop for all particles in the bin (from the last element)
				for (size_t m=0; m<Bins[b].size(); ++m)
				{
					size_t n = Bins[b].size()-1-m;
					size_t p = Bins[b][n];
					if (p>5)
					{
						size_t i = ceil(Lp[p]->X(0)/Lx(0))+1;
						size_t j = ceil(Lp[p]->X(1)/Lx(1))+1;
						size_t k = ceil(Lp[p]->X(2)/Lx(2))+1;

						Lp[p]->CID = i+j*Ncy+k*Ncz;			// bin id
						if (Lp[p]->CID != b)
						{
							swap(Bins[b][n],Bins[b].back());
							Bins[b].pop_back();					// remove particle from bin
							size_t tid = omp_get_thread_num();	// thread id
							updateLpProc[tid].push_back(p);
						}
					}
				}
			}
		}
	}
	
	// for (size_t b=0; b<Bins.size(); ++b)
	// {
	// 	if (Bins[b].size()>0)
	// 	{
	// 		bool stop = false;
	// 		// Loop for all particles in the bin (from the last element)
	// 		for (size_t m=0; m<Bins[b].size(); ++m)
	// 		{
	// 			size_t n = Bins[b].size()-1-m;
	// 			size_t p = Bins[b][n];
	// 			if (p>5)
	// 			{
	// 				size_t i = ceil(Lp[p]->X(0)/Lx(0))+1;
	// 				size_t j = ceil(Lp[p]->X(1)/Lx(1))+1;
	// 				size_t k = ceil(Lp[p]->X(2)/Lx(2))+1;

	// 				Lp[p]->CID = i+j*Ncy+k*Ncz;			// bin id
	// 				if (Lp[p]->CID != b)
	// 				{
	// 					swap(Bins[b][n],Bins[b].back());
	// 					Bins[b].pop_back();					// remove particle from bin
	// 					size_t tid = omp_get_thread_num();	// thread id
	// 					updateLpProc[tid].push_back(p);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// Find the bin list that include less particles
	size_t minP = Lp.size();
	size_t minC = 0;
	for (size_t c=0; c<Nproc; ++c)
	{
		if (BinNump[c]<minP)
		{
			minP = BinNump[c];
			minC = c;
		}
	}

	// Add particles to their bin
	for (size_t c=0; c<Nproc; ++c)
	for (size_t i=0; i<updateLpProc[c].size(); ++i)
	{
		size_t p = updateLpProc[c][i];
		size_t b = Lp[p]->CID;
		Bins[b].push_back(p);
		// Add new actived bin to the bin list that have minmum particles (to balance loading)
		if (Bins[b].size()==1)
		{
			BinsProc[minC].push_back(b);
			// cout << "add bin" << endl;
			// abort(); 
		}
	}

	for (size_t p=6; p<Lp.size(); ++p)
	{
		bool found = false;
		size_t b = Lp[p]->CID;
		for (size_t i=0; i<Bins[b].size(); ++i)
		{
			if (Bins[b][i]==p)
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			cout << "cannot find " << endl;
			abort();
		}
	}		

	vector<vector<vector<size_t>>> lcProc(Nproc);
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=6; p<Lp.size(); ++p)
	{
		size_t tid = omp_get_thread_num();
		bool walls[6] {false, false, false, false, false, false};
		for (size_t i=0; i<NeiV.size(); ++i)
		{
			int b = Lp[p]->CID+NeiV[i];			// neighor bins
			for (size_t j=0; j<Bins[b].size(); ++j)
			{
				size_t q = Bins[b][j];
				if (q<6)
				{
					if (!walls[q])
					{
						double delta = Lp[p]->R-(Lp[p]->X-Lp[q]->X).dot(Lp[q]->Normal);
						if (delta>0.)
						{
							lcProc[tid].push_back({q,p});
							walls[q] = true;
						}
					}
				}
				else if (q>p)
				{
					double delta = Lp[p]->R + Lp[q]->R - (Lp[p]->X-Lp[q]->X).norm();
					if (delta>0.)	lcProc[tid].push_back({p,q});
				}
			}
		}
	}
	// cout << "44444" << endl;
	Lc.clear();
	for (size_t c=0; c<Nproc; ++c)	Lc.insert(Lc.end(), lcProc[c].begin(),lcProc[c].end());
		// cout << "555555555" << endl;
}

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
    	size_t nproc = Nproc;
    	vector<bool> contacts(Lc.size());
    	vector<size_t> keys(Lc.size());
    	vector<Vector3d> xis(Lc.size());
    	vector<Vector3d> xirs(Lc.size());
    	vector<Vector3d> finalPs(Lc.size());

    	#pragma omp parallel for schedule(static) num_threads(nproc)
		for (size_t l=0; l<Lc.size(); ++l)
		{
			int i = Lc[l][0];
			int j = Lc[l][1];

			bool contacted = false;
			Vector3d xi (0.,0.,0.);
			Vector3d xir (0.,0.,0.);
			Vector3d finalP (0.,0.,0.);
			Contact2P(Lp[i], Lp[j], xi, xir, finalP, contacted);

			contacts[l] = contacted;
			keys[l] 	= Key(i,j);
			xis[l]		= xi;
			xirs[l]		= xir;
			finalPs[l] 	= finalP;

			if (Lp[i]->Type==4 || Lp[i]->Type==j)
			for (size_t p=0; p<Lp[i]->MetaP.size(); ++p)
			{
				double ci = CalMetaC(Lp[j]->MetaP, Lp[j]->MetaK, Lp[i]->MetaP[p]);
				if (ci>1.)
				{
					cout << "very big Overlapping between metaballs" << endl;
					cout << "i= " <<i << " j= " << j << endl;
					cout << "tag= " <<Lp[i]->Tag << " tag= " << Lp[j]->Tag << endl;
					// abort();
				}
			}
		}

		FMap.clear();
		RMap.clear();
		MMap.clear();

		for (size_t l=0; l<Lc.size(); ++l)
		{
			if (contacts[l])
			{
				FMap.insert({keys[l], xis[l]});
				RMap.insert({keys[l], xirs[l]});
			}
			MMap.insert({keys[l], finalPs[l]});
		}
		if (writeFc)	WriteContactForceFileH5(n);
    }
}

// inline void DEM::Solve(int tt, int ts, double dt, bool writefile)
// {
// 	double time0, time1, time2;
// 			time0 = 0.;
// 			time1 = 0.;
// 			time2 = 0.;
// 	Dt = dt;
// 	for (int t=0; t<tt; ++t)
// 	{
// 		bool show = false;
// 		// bool show = true;
// 		// if (t%ts==0)	show = true;
// 		bool save = false;
// 		if (t%ts == 0)
// 		{
// 			save = true;
// 			/*if (t>250000)*/cout << "Time Step ============ " << t << endl;
// 			// if (writefile /*&& t>4000*/)	WriteFileH5(t);
// 			// time0 = 0.;
// 			// time1 = 0.;
// 			// time2 = 0.;
// 		}
// 		if (t>347000 && t<348000)
// 		{
// 			cout << "Time Step ============ " << t << endl;
// 			cout << "x: " << Lp[6]->X.transpose() << endl;
// 			cout << "v: " << Lp[6]->V.transpose() << endl;
// 		}
// 		// if (t>405500)
// 		// {		size_t i,j,k;
// 		// 		FindIndex(Lp[151]->CID, i, j, k);
// 		// 		cout << "p151: " << Lp[151]->CID << " i: " << i << " j: " <<j << " k: " << k << endl;
// 		// 		FindIndex(Lp[155]->CID, i, j, k);
// 		// 		cout << "p155: " << Lp[155]->CID << " i: " << i << " j: " <<j << " k: " << k << endl;
// 		// }		// double lm = 0.5*Lp[6]->M*Lp[6]->V.squaredNorm() + 0.5*Lp[7]->M*Lp[7]->V.squaredNorm();
// /*		Vector3d lmi, lmj, lm;
// 		lmi = Lp[6]->M*Lp[6]->V;
// 		lmj = Lp[7]->M*Lp[7]->V;
// 		lm = lmi+lmj;*/
// 		// cout << "Lp[6]->M: " << Lp[6]->M << endl;
// 		// cout << "v: " << Lp[6]->V.transpose() << endl;
// 		// cout << "Lp[7]->M: " << Lp[7]->M << endl;
// 		// cout << "v: " << Lp[7]->V.transpose() << endl;
// 		Vector3d vi0 (0., 0.001, 0.);
// 		Vector3d vj0 (0., -0.001, 0.);
// 		Vector3d lmi0 = Lp[6]->M*vi0;
// 		Vector3d lmj0 = Lp[7]->M*vj0;
// 		cout << "lniear momentum at beginning: " << (lmi0+lmj0).transpose() << endl;
// 		cout << "lniear momentum: " << lm.transpose() << endl;
// 		// cout << (Lp[6]->M*Lp[6]->V + Lp[7]->M*Lp[7]->V).transpose() << endl;
// 		// abort();
// 		Vector3d ami, amj, am;
// 		ami = Lp[6]->X.cross(lmi);
// 		amj = Lp[7]->X.cross(lmj);
// 		am = ami + amj;
// 		// ami << Lp[6]->I(0)*Lp[6]->W(0), Lp[6]->I(1)*Lp[6]->W(1), Lp[6]->I(2)*Lp[6]->W(2);
// 		// amj << Lp[7]->I(0)*Lp[7]->W(0), Lp[7]->I(1)*Lp[7]->W(1), Lp[7]->I(2)*Lp[7]->W(2);
// 		// am = Lp[6]->Q._transformVector(ami) + Lp[7]->Q._transformVector(amj);

// 		cout << "angular momentum at beginning: " << 0. << endl;
// 		cout << "angular momentum: " << am.transpose() << endl;
// 		// abort();

// 		double en = Lp[6]->CalEnergy() + Lp[7]->CalEnergy();
// 		cout << "Energy at beginning: " << 0.5*0.001*0.001*(Lp[6]->M + Lp[7]->M) << endl;
// 		cout << "Energy: " << en << endl;
// 		// abort();
// 		// if (show)	cout << "Time Step ============ " << t << endl;
// 		auto t_start = std::chrono::system_clock::now();
// 		// FindContact();
// 		// cout << "LinkedCell start" << endl;
// 		LinkedCell();
// 		// cout << "LinkedCell end" << endl;
// 		// Lc.push_back({6,7});
// 		// bool firstStep = false;
// 		// if (t==0)	firstStep = true;
// 		// FindContactBasedOnNode(firstStep);
// 		auto t_end = std::chrono::system_clock::now();
// 		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		// if (show)	cout << "FindContact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 		if (show)
// 		{
// 			cout << "FindContact time= " << time0 << endl;
// 			time0 = 0.;
// 		}
// 		t_start = std::chrono::system_clock::now();
// 		// cout << "start contact" << endl;
// 		Contact(false, 0);
// 		// cout << "Finish contact" << endl;
// 		int nc = Lc.size();
// 		Lc.clear();
// 		t_end = std::chrono::system_clock::now();
// 		if (show)   cout << "nc= " << nc << endl;
// 		time1 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		// if (show)	cout << "Contact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 		if (show)
// 		{
// 			cout << "Contact time= " << time1 << endl;
// 			time1 = 0.;
// 		}
// 		if (t==0)
// 		{
// 			for (size_t p=0; p<Lp.size(); ++p)
// 			{
// 				Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
// 				Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
// 			}
// 		}
// 		t_start = std::chrono::system_clock::now();

// 		Move(save);
// 		// bool movep = false;
// 		// if (t>405500 && t%100==0)	movep = true;
// 		// Move(movep);
// 		t_end = std::chrono::system_clock::now();
// 		time2 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
// 		// if (show)	cout << "Move time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 		if (show)
// 		{
// 			cout << "Move time= " << time2 << endl;
// 			time2 = 0.;
// 		}
// 		t_start = std::chrono::system_clock::now();
// 		// WriteFileParticleInfo(t);
// 		ZeroForceTorque(true, true);
// 		t_end = std::chrono::system_clock::now();
// 		// if (show)	cout << "ZeroForceTorque time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
// 		// if (writefile && save /*&& t>111000*/)	WriteFileH5(t);
// 		if (t>347000 && t<348000)	WriteFileH5(t);
// 	}
// }

inline void DEM::SolveOneStep(double dt, bool save)
{
	UpdatePositionGlobal(dt, save, true);
	LinkedCell();
	Contact(false, 0);
	UpdateVelocityGlobal(dt);
	Lc.clear();
}

inline void DEM::Solve(int tt, int ts, double dt, bool writefile)
{
	double time0, time1, time2;
			time0 = 0.;
			time1 = 0.;
			time2 = 0.;
	Dt = dt;
	GenBinsProc();
	for (int t=0; t<tt; ++t)
	{
		bool show = false;
		show = true;
		// if (t%ts==0)	show = true;
		bool save = false;
		if (t%ts == 0)
		{
			save = true;
			/*if (t>250000)*/cout << "Time Step ============ " << t << endl;
			// if (writefile /*&& t>4000*/)	WriteFileH5(t);
			// time0 = 0.;
			// time1 = 0.;
			// time2 = 0.;
		}
		// if (t>347000 && t<348000)
		// {
		// 	cout << "Time Step ============ " << t << endl;
		// 	cout << "x: " << Lp[6]->X.transpose() << endl;
		// 	cout << "v: " << Lp[6]->V.transpose() << endl;
		// }
		if (t==0)
		{
			for (size_t p=0; p<Lp.size(); ++p)
			{
				Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
				Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
			}
		}
		// if (show)	cout << "Time Step ============ " << t << endl;
		auto t_start = std::chrono::system_clock::now();
		UpdatePositionGlobal(dt, save, true);
		auto t_end = std::chrono::system_clock::now();
		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
		if (show)
		{
			cout << "UpdatePositionGlobal time= " << time0 << endl;
			time0 = 0.;
		}

		t_start = std::chrono::system_clock::now();
		LinkedCell();
		t_end = std::chrono::system_clock::now();
		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
		if (show)
		{
			cout << "LinkedCell time= " << time0 << endl;
			time0 = 0.;
		}

		t_start = std::chrono::system_clock::now();
		Contact(false, 0);
		t_end = std::chrono::system_clock::now();
		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
		if (show)
		{
			cout << "Contact time= " << time0 << endl;
			time0 = 0.;
		}

		t_start = std::chrono::system_clock::now();
		UpdateVelocityGlobal(dt);
		t_end = std::chrono::system_clock::now();
		time0 += std::chrono::duration<double, std::milli>(t_end-t_start).count();
		if (show)
		{
			cout << "UpdateVelocityGlobal time= " << time0 << endl;
			time0 = 0.;
		}
		int nc = Lc.size();
		Lc.clear();
		if (show)   cout << "nc= " << nc << endl;
		if (writefile && save /*&& t>111000*/)	WriteFileH5(t);
		// if (t>83000 && t<84000)	WriteFileH5(t);
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

inline void DEM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(9) << setfill('0') << n;			
	string file_name_h5 = "DEM_"+out.str()+".h5";

    hid_t     file_id;
    file_id = H5Fcreate(file_name_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);	
	
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

	float* r_h5 	= new float[  npar];
	float* rho_h5 	= new float[  npar];
	float* tag_h5 	= new float[  npar];
	float* pos_h5 	= new float[3*npar];
	float* vel_h5 	= new float[3*npar];
	float* agv_h5	= new float[3*npar];
	float* fh_h5 	= new float[3*npar];
	float* poi_h5	= new float[3*n_points];
	int* 	fac_h5	= new int   [n_fe];
	float* fv_h5 	= new float[n_faces];
	float* ftag_h5 = new float[n_faces];
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

	H5LTmake_dataset_float(file_id,"Radius",1,dims_scalar,r_h5);
	H5LTmake_dataset_float(file_id,"Rho",1,dims_scalar,rho_h5);
	H5LTmake_dataset_float(file_id,"Tag",1,dims_scalar,tag_h5);
	H5LTmake_dataset_float(file_id,"Position",1,dims_vector,pos_h5);
	H5LTmake_dataset_float(file_id,"Velocity",1,dims_vector,vel_h5);
	H5LTmake_dataset_float(file_id,"AngularVelocity",1,dims_vector,vel_h5);
	H5LTmake_dataset_float(file_id,"HydroForce",1,dims_vector,fh_h5);
	H5LTmake_dataset_float(file_id,"Points",1,dims_points,poi_h5);
	H5LTmake_dataset_int(file_id,"Faces",1,dims_faces,fac_h5);
	H5LTmake_dataset_float(file_id,"FaceVelocity",1,dims_fscalar,fv_h5);
	H5LTmake_dataset_float(file_id,"FaceTag",1,dims_fscalar,ftag_h5);

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

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

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