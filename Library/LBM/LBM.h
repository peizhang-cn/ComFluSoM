#include "../HEADER.h"
#include <cstring>

// discrete model
enum DnQm
{
	D2Q9 ,
	D3Q15,
    D3Q19,
    D3Q27
};

enum CollisionModel
{
	SRT,
	MRT,
    CM
};

class LBM
{
public:
	LBM();
	~LBM();
	LBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu);
	void Init(double rho0, Vector3d initV);
	void CalRhoVLocal(int i, int j, int k);													// Calculate density and velocity for local lattice
	void CalRhoV();																			// Calculate density and velocity
	void CalRho();																			// Calculate only density
	void Stream();
	void BodyForceLocal(int i, int j, int k, Vector3d force);								// Calculate body force
	void BodyForceLocalOpenMP(int i, int j, int k, Vector3d force);							// Calculate body force OpenMP version
	void BodyForceLocalTwoStep(int i, int j, int k, Vector3d force);
	void BounceBack(int i, int j, int k);
	void SBounceBack(int i, int j, int k);
	void Boundary(double delta, int bctype, Vector3d vw);
	void SetA(Vector3d a);
	void SetPeriodic(bool x, bool y, bool z);
	// void WriteFileH5(int n);
	void WriteFileH5(int n, int scale);
	void BoundaryAC(Vector3d vb);
	void VIBB(int i, int j, int k, int q, double delta, Vector3d& vw, Vector3d& fh);
	void IBBYu(int i, int j, int k, int q, double delta, Vector3d& vw, Vector3d& fh);
	void VAM(int i, int j, int k, double g, double rhos, Vector3d& vs, Vector3d& fh);
	void PSM(int i, int j, int k, double g, Vector3d& vw, Vector3d& fh);
	void ApplyWall();
	void SetWall();
	void FindIndex(int n, int& i, int& j, int& k);
	void ReadG(string fileName);
	Vector3d InterpolateV(const Vector3d& x);
	void ActiveLES(double sc);
	/*===================================Methods for SRT=====================================================*/
	void  (LBM::*CalFeq)(VectorXd& feq, double phi, Vector3d v);							// Function pointer to calculate equilibrium distribution 
	void CalFeqC(VectorXd& feq, double phi, Vector3d v);									// Function to calculate equilibrium distribution for compressible fluid
	void CalFeqIC(VectorXd& feq, double phi, Vector3d v);									// Function to calculate equilibrium distribution for incompressible fluid
	void CollideSRTLocal(int i, int j, int k);
	void CollideSRTLocal(int i, int j, int k, double nu);									// For variable viscosity
	void CollideSRT();
	/*===================================Methods for MRT=====================================================*/
	void  (LBM::*CalMeq)(VectorXd& meq, double rho, Vector3d v);							// Function pointer to calculate equilibrium distribution 
	void CalMeqD2Q9(VectorXd& meq, double rho, Vector3d v);									// Function to calculate equilibrium distribution for compressible fluid
	void CalMeqD3Q15(VectorXd& meq, double rho, Vector3d v);
	void CollideMRTLocal(int i, int j, int k);
	void CollideMRT();
	/*===================================Methods for CM =====================================================*/
	int 							Nproc;													// Number of processors which used

	double*** 						Rho;													// Fluid density
	double***						Omega;													// Reversed Tua for variable viscosities
	double****	 					G;														// Flag of lattice type
	vector<size_t>***		 		Flag;													// Flag of lattice type
	Vector3d***						V;														// Fluid velocity
	Vector3d***						ExForce;												// External force

	Vector3d						A;														// Global Acceleration

// private:

//Private Variables =====================================================================================================================================

	CollisionModel 					Cmodel;
    int 							Nx;														// Domain size
    int 							Ny;
    int 							Nz;
	int 							DomSize[3];    
    int 							Ncell;													// Number of lattices
    int 							Ncz;
    int 							Ncy;

    double 							Nu;														// Vsicosity
    double 							Tau;													// Relaxation time                                                     
    double 							Omega0;													// Reciprocal of relaxation time
    double 							Rho0;
    double 							Sc2;													// Square of Smagorinsky constant

    int 							D;														// Dimension
    int 							Q;														// Number of discrete velocity vectors
    vector<Vector3d> 				E;														// Discrete velocities
    vector<Vector3i> 				Ne;														// Relative location of neighbor cells
    vector <double> 				W;														// Weights
    vector <int> 					Op;														// Opposite directions

    MatrixXd						S;														// Relaxation matrix
    MatrixXd						M;														// Transform matrix
    MatrixXd						Mi;														// Inverse of transform matrix
    MatrixXd						Ms;														// Inverse of M multiply by S, for speed up.
    MatrixXd						Mf;														// Inverse of M multiply by (I-0.5*S) for force term, for speed up

    VectorXd*** 					F;														// Distribution function
    VectorXd*** 					Ft;														// Distribution function

    vector<Vector3i>				Lwall;													// List of wall nodes
    bool 							InCompressible;											// Whether the fluid is incompressible or not 
    bool 							Periodic[3];											// Whether periodic on x, y and z direction
    bool 							Convergence;											// Whether convergence for steady problems
};

inline LBM::LBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu)
{
	Cmodel = cmodel;
	Nx	= nx;
	Ny	= ny;
	Nz	= nz;
	DomSize[0] = Nx;
	DomSize[1] = Ny;
	DomSize[2] = Nz;

	Ncell = (Nx+1)*(Ny+1)*(Nz+1);

	Ncz = (Nx+1)*(Ny+1);
	Ncy = (Nx+1);

	Nu  = nu;
    Tau = 3.*Nu + 0.5;
    Omega0 = 1./Tau;
    Sc2 = 0.;
    cout << "Relaxation time = " << Tau << endl;
    Nproc = 12;

	if (dnqm==D2Q9)
	{
		E   = {     { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0}, 
                	{-1, 0, 0}, { 0,-1, 0}, { 1, 1, 0}, 
                	{-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
		W   = {     4./9. , 1./9. , 1./9. ,
                	1./9. , 1./9. , 1./36., 
                	1./36., 1./36., 1./36. };
		Op  = {     0, 3, 4, 
                	1, 2, 7, 
                	8, 5, 6 };
		Ne  = {     { 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };
		D   = 2;
		Q   = 9;

	    if (Cmodel==MRT)
	    {
		    M.resize(Q,Q);
		    S.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0) <<  1,  1,  1,  1,  1,  1,  1,  1,  1;
		    M.row(1) << -4, -1, -1, -1, -1,  2,  2,  2,  2;
		    M.row(2) <<  4, -2, -2, -2, -2,  1,  1,  1,  1;
		    M.row(3) <<  0,  1,  0, -1,  0,  1, -1, -1,  1;
		    M.row(4) <<  0, -2,  0,  2,  0,  1, -1, -1,  1;
		    M.row(5) <<  0,  0,  1,  0, -1,  1,  1, -1, -1;
		    M.row(6) <<  0,  0, -2,  0,  2,  1,  1, -1, -1;
		    M.row(7) <<  0,  1, -1,  1, -1,  0,  0,  0,  0;
		    M.row(8) <<  0,  0,  0,  0,  0,  1, -1,  1, -1;

		    S = MatrixXd::Identity(Q,Q);
		    S(1,1) = 0.3;
		    S(2,2) = 1.5;
		    S(4,4) = S(6,6) = 1.2;
		    S(7,7) = S(8,8) = Omega0;

		    Mi = M.inverse();
		    Ms = M.inverse()*S;
		    Mf = M.inverse()*(MatrixXd::Identity(Q,Q) - 0.5*S);

		    CalMeq = &LBM::CalMeqD2Q9;
	    }
	}

	else if (dnqm==D3Q15)
	{
		E   = {	   { 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
			       { 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
			       {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, {-1, 1, 1}, { 1,-1,-1} };
		W   = {	   2./9. , 1./9. , 1./9. , 1./9. , 1./9. ,
			       1./9. , 1./9. , 1./72., 1./72., 1./72.,
			       1./72., 1./72., 1./72., 1./72., 1./72. };
		Op  = {	   0,  2,  1,  4,  3, 
                   6,  5,  8,  7, 10, 
                   9, 12, 11, 14, 13 };
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
		D   = 3;
		Q   = 15;

	    if (Cmodel==MRT)
	    {
		    M.resize(Q,Q);
		    S.resize(Q,Q);
		    Ms.resize(Q,Q);
		    Mf.resize(Q,Q);

		    M.row(0 ) <<  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(1 ) << -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(2 ) << 16.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0;
		    M.row(3 ) <<  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0;
		    M.row(4 ) <<  0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0;
		    M.row(5 ) <<  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0;
		    M.row(6 ) <<  0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0;
		    M.row(7 ) <<  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0;
		    M.row(8 ) <<  0.0,  0.0,  0.0,  0.0,  0.0, -4.0,  4.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0;
		    M.row(9 ) <<  0.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0;
		    M.row(10) <<  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0;
		    M.row(11) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0;
		    M.row(12) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0;
		    M.row(13) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0;
		    M.row(14) <<  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0;

		    S = MatrixXd::Identity(Q,Q);
		    S(0,0) = 0.;
		    S(1,1) = 1.6;
		    S(2,2) = 1.2;
		    S(3,3) = S(5,5) = S(7,7) = 0.;
		    S(4,4) = S(6,6) = S(8,8) = 1.6;
		    S(5,5) = 0.;
		    S(7,7) = 0.;
		    S(9,9) = S(10,10) = S(11,11) = S(12,12) = S(13,13) = Omega0;
		    S(14,14) = 1.2;

		    Mi = M.inverse();
		    Ms = M.inverse()*S;
		    Mf = M.inverse()*(MatrixXd::Identity(Q,Q) - 0.5*S);

		    CalMeq = &LBM::CalMeqD3Q15;
	    }
	}

    else if (dnqm==D3Q19)
    {
        E   = {    { 0, 0, 0}, 
                   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
                   { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},
                   { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1} };
        W   = {    1./3. ,
                   1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 
                   1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 
                   1./36., 1./36., 1./36., 1./36., 1./36., 1./36. };
        Op  = {    0 , 
                   2 , 1 , 4 , 3 , 6 , 5 , 
                   8 , 7 , 10, 9 , 12, 11, 
                   14, 13, 16, 15, 18, 17 };
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
        D   = 3;
        Q   = 19;
    }
    else if (dnqm==D3Q27)
    {
    	// // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    	// // 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1
    	//    0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1
    	// // 0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1
    	//    0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1
    	// // 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1
    	//    0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1

        E   = {    { 0, 0, 0},

                   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},

                   { 1, 1, 0}, {-1, 1, 0}, { 1,-1, 0}, {-1,-1, 0}, 
                   { 1, 0, 1}, {-1, 0, 1}, { 1, 0,-1}, {-1, 0,-1},
                   { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1},

                   { 1, 1, 1}, {-1, 1, 1}, { 1,-1, 1}, {-1,-1, 1},
                   { 1, 1,-1}, {-1, 1,-1}, { 1,-1,-1}, {-1,-1,-1} };
        W   = {    8./27., 
        		   2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 
        		   1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 
        		   1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.};
        Op  = {    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19};
		Ne  = {	   { 0, 0, 0},
				   { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
		           { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}, 
		           { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
			       { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}, 
			       { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} };
        D   = 3;
        Q   = 27;
    }
    else
    {
        cout << "Undefined DnQm model: " << dnqm << endl;
        abort();
    }
	
	InCompressible = incompressible;
	if (InCompressible)
	{
		CalFeq = &LBM::CalFeqIC;
		cout << "Using incompressible Feq, maybe not a good idea" << endl;
	}
	else					CalFeq = &LBM::CalFeqC;

	Convergence = false;
	Periodic[0] = true;
	Periodic[1] = true;
	Periodic[2] = true;

	A = Vector3d::Zero();
	Lwall.resize(0);
}

inline void LBM::Init(double rho0, Vector3d initV)
{
    cout << "================ Start init. ================" << endl;

    Rho0 	= rho0;
	Rho		= new double**  [Nx+1];
	Omega	= new double**  [Nx+1];
	G		= new double***	[Nx+1];
	Flag 	= new vector<size_t>** [Nx+1];
	V		= new Vector3d**[Nx+1];
	ExForce	= new Vector3d**[Nx+1];
	F		= new VectorXd**[Nx+1];
	Ft		= new VectorXd**[Nx+1];

	for (int i=0; i<=Nx; ++i)
	{
		Rho    [i]	= new double*  [Ny+1];
		Omega  [i]	= new double*  [Ny+1];
		G      [i]	= new double** [Ny+1];
		Flag   [i]  = new vector<size_t>* [Ny+1];
		V      [i]	= new Vector3d*[Ny+1];
		ExForce[i]	= new Vector3d*[Ny+1];
		F      [i]	= new VectorXd*[Ny+1];
		Ft     [i]	= new VectorXd*[Ny+1];

		for (int j=0; j<=Ny; j++)
		{
			Rho    [i][j]	= new double  [Nz+1];
			Omega  [i][j]	= new double  [Nz+1];
			G      [i][j]	= new double* [Nz+1];
			Flag   [i][j]	= new vector<size_t> [Nz+1];
			V      [i][j]	= new Vector3d[Nz+1];
			ExForce[i][j]	= new Vector3d[Nz+1];
			F      [i][j]	= new VectorXd[Nz+1];
			Ft     [i][j]	= new VectorXd[Nz+1];

			for (int k=0; k<=Nz; k++)
			{
				Rho    [i][j][k]	= Rho0;
				Omega  [i][j][k]	= Omega0;
				G      [i][j][k]	= new double[4];
				Flag   [i][j][k].resize(0);
				V      [i][j][k]	= initV;
				// V      [i][j][k]	= 4.*initV*((j-1)/(double) (Ny-2))*(1.-(j-1)/(double) (Ny-2));
				ExForce[i][j][k]	= Vector3d::Zero();

				G[i][j][k][0] = -1.;
				G[i][j][k][1] = 0.;
				G[i][j][k][2] = -1.;
				G[i][j][k][3] = 0.;

				VectorXd feq(Q);
				if (Cmodel==SRT)
				{
					(this->*CalFeq)(feq, rho0, initV);
				}
				else if (Cmodel==MRT)
				{
					VectorXd meq(Q);
					(this->*CalMeq)(meq, rho0, initV);
					feq = Mi*meq;
				}

				F [i][j][k]	= feq;
				Ft[i][j][k]	= feq;
			}
		}
	}

	cout << "================ Finish init. ================" << endl;
}

LBM::ActiveLES(double sc)
{
	if (sc<0.1 || sc>0.2)
	{
		cout << "Smagorinsky constant should be within range of 0.1~0.2." << endl;
		abort();
	}
	Sc2 = sc*sc;
}

// inline void LBM::ReadG(string fileName)
// {
// 	cout << fileName << endl;
// 	ifstream file(fileName);
// 	// cout << "open file" << endl;

// 	if(file.is_open())
// 	{
//         int i,j,k,g;
//         while (file >> i >> j >> k >> g){
//         	if (g==1)		G[i][j][k] = -2.;
//         }
//         // for (int i = 0; i <= Nx; ++i)
//         // for (int j = 0; j <= Ny; ++j)
//         // for (int k = 0; k <= Nz; ++k)
//         // {
//         // 	// cout << "i= " << i << " j= " << j << " k= " << k << endl;
//         // 	file >> G[i][j][k];

//         // 	if (G[i][j][k]==1.)		G[i][j][k] = -2.;

//         // 	// cout << "g= " << G[i][j][k] << endl;
//         // }	
// 	}
// 	else
// 	{
// 		cout << "Cannot find " << fileName << "!!!" << endl;
// 		abort();
// 	}

// 	cout << "Finish ReadG" << endl;

//     // for (int i = 0; i <= Nx; ++i)
//     // for (int j = 0; j <= Ny; ++j)
//     // for (int k = 0; k <= Nz; ++k)
//     // {
//     // 	for (size_t q=0; q<Ne.size(); ++q)
//     // 	{
//     // 		int in = i + Ne[q][0];
//     //         int jn = j + Ne[q][1];
//     //         int kn = k + Ne[q][2];
//     // 	}
//     // }
// }

inline void LBM::CalFeqC(VectorXd& feq, double rho, Vector3d v)
{
    double vv  = v.dot(v);
    for (int q=0; q<Q; ++q)
    {
    	double ev  = E[q].dot(v);
    	feq(q) = W[q]*rho*(1. + 3.*ev + 4.5*ev*ev - 1.5*vv);
	}
}

inline void LBM::CalFeqIC(VectorXd& feq, double rho, Vector3d v)
{
    double vv  = v.dot(v);
    double deltaRho = rho-Rho0;
    for (int q=0; q<Q; ++q)
    {
    	double ev  = E[q].dot(v);
    	feq(q) = W[q]*(deltaRho + Rho0*(1. + 3.*ev + 4.5*ev*ev - 1.5*vv));
	}
}

inline void LBM::CalMeqD2Q9(VectorXd& meq, double rho, Vector3d v)
{
	meq(0) = rho;
	meq(1) = rho*(-2 + 3*(v(0)*v(0) + v(1)*v(1)));
	meq(2) = -rho - meq(1);
	meq(3) = rho*v(0);
	meq(4) = -meq(3);
	meq(5) = rho*v(1);
	meq(6) = -meq(5);
	meq(7) = rho*(v(0)*v(0) - v(1)*v(1));
	meq(8) = rho*v(0)*v(1);
}

inline void LBM::CalMeqD3Q15(VectorXd& meq, double rho, Vector3d v)
{
	double v2 = v.squaredNorm();
	meq(0) = rho;
	meq(1) = rho*(v2-1.);
	meq(2) = rho*(1.-5*v2);
	meq(3) = rho*v(0);
	meq(4) = -7.*rho*v(0)/3.;
	meq(5) = rho*v(1);
	meq(6) = -7.*rho*v(1)/3.;
	meq(7) = rho*v(2);
	meq(8) = -7.*rho*v(2)/3.;
	meq(9) = rho*(3.*v(0)*v(0)-v2);
	meq(10)= rho*(v(1)*v(1)-v(2)*v(2));
	meq(11)= rho*v(0)*v(1);
	meq(12)= rho*v(1)*v(2);
	meq(13)= rho*v(0)*v(2);
	meq(14)= 0.;
}

inline void LBM::CalRhoVLocal(int i, int j, int k)
{
	Rho[i][j][k] = 0.;
    V[i][j][k] = Vector3d::Zero();

    for (int q = 0; q < Q; ++q)
    {
    	Rho[i][j][k]	+= F[i][j][k](q);
    	V[i][j][k] 		+= F[i][j][k](q)*E[q];
    }

    V[i][j][k] /= Rho[i][j][k];
}

inline void LBM::CalRhoV()
{
	// double err = 0.;
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i = 0; i <= Nx; ++i)
    for (int j = 0; j <= Ny; ++j)
    for (int k = 0; k <= Nz; ++k)
    {
		CalRhoVLocal(i,j,k);
    }
}

inline void LBM::CalRho()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i = 0; i <= Nx; ++i)
    for (int j = 0; j <= Ny; ++j)
    for (int k = 0; k <= Nz; ++k)
    {
    	Rho[i][j][k] = F[i][j][k].sum();
    }	
}

inline void LBM::CollideSRTLocal(int i, int j, int k)
{
	VectorXd feq(Q), fneq(Q);
	(this->*CalFeq)(feq, Rho[i][j][k], V[i][j][k]);

	fneq = F[i][j][k] - feq;

	double OmegaEddy = 0.;
	if (Sc2>0.)
	{
		double qHat2 = 0.;
		for (size_t m=0; m<D; ++m)
		for (size_t n=0; n<D; ++n)
		{
			double qij = 0.;
			for (size_t q=0; q<Q; ++q)	qij += E[q](m)*E[q](n)*fneq(q);
			qHat2 += qij*qij;
		}
		double qHat = 1.41414141414*sqrt(qHat2);
		OmegaEddy = 2./(sqrt(Tau+18.*Sc2*qHat)-Tau);
	}

	Ft[i][j][k] = F[i][j][k] + (Omega[i][j][k]+OmegaEddy)*Fneq;
}

inline void LBM::CollideMRTLocal(int i, int j, int k)
{
	VectorXd meq(Q), m(Q);
	(this->*CalMeq)(meq, Rho[i][j][k], V[i][j][k]);
	m = M*F[i][j][k];

    Ft[i][j][k] = F[i][j][k] + Ms*(meq - m);
}

// Body force (Eq.14) 
// Unified Theory of Lattice Boltzmann Models for Nonideal Gases
inline void LBM::BodyForceLocal(int i, int j, int k, Vector3d force)
{
	for (int q=0; q<Q; ++q)
	{
		Vector3d eeu = E[q] - V[i][j][k] + 3*E[q].dot(V[i][j][k])*E[q];
		Ft[i][j][k](q) += 3*W[q]*eeu.dot(force);

		if (Ft[i][j][k](q)<0.)
		{
			cout << "negative after force term!" << endl;
			cout << "force= " << force.transpose() << endl;
			cout << i << " " << j << " " << k << endl;
			cout << "vel= " << V[i][j][k].transpose() << endl;
			abort();
		}
	}
}

inline void LBM::BodyForceLocalOpenMP(int i, int j, int k, Vector3d force)
{
	for (int q=0; q<Q; ++q)
	{
		Vector3d eeu = E[q] - V[i][j][k] + 3*E[q].dot(V[i][j][k])*E[q];
		double fb = 3*W[q]*eeu.dot(force);
		#pragma omp atomic
		Ft[i][j][k](q) += fb;

		if (Ft[i][j][k](q)<0.)
		{
			cout << "negative after force term!" << endl;
			cout << "force= " << force.transpose() << endl;
			cout << i << " " << j << " " << k << endl;
			cout << "vel= " << V[i][j][k].transpose() << endl;
			abort();
		}
	}
}

inline void LBM::BodyForceLocalTwoStep(int i, int j, int k, Vector3d force)
{
	VectorXd fb(Q);

	for (int q=0; q<Q; ++q)
	{
		Vector3d eeu = 3.*(E[q] - V[i][j][k]) + 9.*E[q].dot(V[i][j][k])*E[q];
		fb(q) = W[q]*(1-0.5*Omega[i][j][k])*eeu.dot(force);
	}

	Ft[i][j][k] += fb;
	for (int q=0; q<Q; ++q)
	{
		if (Ft[i][j][k](q)<0.)
		{
			cout << "negative after force term!" << endl;
			cout << "force= " << force.transpose() << endl;
			cout << i << " " << j << " " << k << endl;
			cout << "vel= " << V[i][j][k].transpose() << endl;
			abort();
		}
	}
}


// inline void LBM::CollideSRT()
// {
//     // cout << "CollideSRT" << endl;
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (int i=0; i<=Nx; i++)
//     for (int j=0; j<=Ny; j++)
//     for (int k=0; k<=Nz; k++)
//     {
//     	CollideSRTLocal(i, j, k);
//     	Vector3d bf = Rho[i][j][k]*A + ExForce[i][j][k];
//     	BodyForceLocal(i, j, k, bf);
//     	ExForce[i][j][k] = Vector3d::Zero();
//     }
// }

inline void LBM::FindIndex(int n, int& i, int& j, int& k)
{
	k = n/Ncz;
	j = (n%Ncz)/Ncy;
	i = (n%Ncz)%Ncy;
}


// inline void LBM::CollideSRT()
// {
//     // cout << "CollideSRT" << endl;
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (int n=0; n<Ncell; ++n)
//     {
//     	int i, j, k;
//     	FindIndex(n, i, j, k);

//     	CollideSRTLocal(i, j, k);
//     	// Vector3d bf = Rho[i][j][k]*A + ExForce[i][j][k];
//     	// BodyForceLocal(i, j, k, bf);
//     	// ExForce[i][j][k] = Vector3d::Zero();
//     }
// }

inline void LBM::CollideSRT()
{
    // cout << "CollideMRT" << endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i=0; i<=Nx; i++)
    for (int j=0; j<=Ny; j++)
    for (int k=0; k<=Nz; k++)
    {
		CollideSRTLocal(i, j, k);
	    Vector3d bf = Rho[i][j][k]*A + ExForce[i][j][k];
    	BodyForceLocal(i, j, k, bf);
	    ExForce[i][j][k] = Vector3d::Zero();
    }
}

inline void LBM::CollideMRT()
{
    // cout << "CollideMRT" << endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int n=0; n<Ncell; ++n)
    {
    	int i, j, k;
    	FindIndex(n, i, j, k);

    	CollideMRTLocal(i, j, k);
    	Vector3d bf = Rho[i][j][k]*A + ExForce[i][j][k];
    	BodyForceLocal(i, j, k, bf);
    	ExForce[i][j][k] = Vector3d::Zero();
    }
}


// inline void LBM::CollideMRT()
// {
//     // cout << "CollideMRT" << endl;
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (int i=0; i<=Nx; i++)
//     for (int j=0; j<=Ny; j++)
//     for (int k=0; k<=Nz; k++)
//     {
//     	if (G[i][j][k]==-2.)	SBounceBack(i,j,k);
// 	    else 					CollideMRTLocal(i, j, k);
// 	    Vector3d bf = Rho[i][j][k]*A + ExForce[i][j][k];
//     	BodyForceLocal(i, j, k, bf);
// 	    ExForce[i][j][k] = Vector3d::Zero();
//     }
// }

inline void LBM::SetA(Vector3d a)
{
    A = a;
}

inline void LBM::SetPeriodic(bool x, bool y, bool z)
{
	Periodic[0] = x;
	Periodic[1] = y;
	Periodic[2] = z;
}

inline void LBM::BounceBack(int i, int j, int k)
{
	for (int q=0; q< Q;  ++q)
	{
		int ip = (i - (int) E[q][0]);
		int jp = (j - (int) E[q][1]);
		int kp = (k - (int) E[q][2]);

		bool out[3] = {false, false, false};

		if (ip<0 || ip>Nx)	out[0] = true;
		if (jp<0 || jp>Ny)	out[1] = true;
		if (kp<0 || kp>Nz)	out[2] = true;

		if (out[0] || out[1] || out[2])
		{
			F[i][j][k](q) 	= Ft[i][j][k](Op[q]);
		}
	}
}

inline void LBM::SBounceBack(int i, int j, int k)
{
	VectorXd ft(Q);
	ft = Ft[i][j][k];
	for (int q=0; q< Q;  ++q)
	{

		Ft[i][j][k](q) 	= ft(Op[q]);
	}
}

inline void LBM::Stream()
{
    // cout << "Stream" << endl;
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	for (int q=0; q< Q;  ++q)
	{
		int ip = (i- (int) E[q][0]+Nx+1)%(Nx+1);
		int jp = (j- (int) E[q][1]+Ny+1)%(Ny+1);
		int kp = (k- (int) E[q][2]+Nz+1)%(Nz+1);

		F[i][j][k][q] = Ft[ip][jp][kp][q];
	}
}

// inline void LBM::Stream()
// {
//     // cout << "Stream" << endl;
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (int n=0; n<Ncell; ++n)
//     {
//     	int i, j, k;
//     	FindIndex(n, i, j, k);
    	
// 		for (int q=0; q< Q;  ++q)
// 		{
// 			int ip = (i- (int) E[q][0]+Nx+1)%(Nx+1);
// 			int jp = (j- (int) E[q][1]+Ny+1)%(Ny+1);
// 			int kp = (k- (int) E[q][2]+Nz+1)%(Nz+1);

// 			F[i][j][k][q] = Ft[ip][jp][kp][q];
// 		}
//     }
// }

inline void LBM::BoundaryAC(Vector3d vb)
{
	// Vector3d v1 (0.001, 0., 0.);
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		// Vector3d v0 = vb;
		// Vector3d v0 = 1.5*vb*j*(Ny-j)/(0.25*Ny*Ny);

		Vector3d v0 = 4.*vb*((j-1)/(double) (Ny-2))*(1.-(j-1)/(double) (Ny-2));

		Rho[0][j][k] = Rho[1][j][k];
		V  [0][j][k] = v0;

		VectorXd fneq(Q), feq(Q);
		(this->*CalFeq)(feq, Rho[1][j][k], V[1][j][k]);
		fneq = F[1][j][k] - feq;
		(this->*CalFeq)(Ft[0][j][k], Rho[0][j][k], V[0][j][k]);
		Ft[0][j][k] += fneq;

		Ft[Nx][j][k] = F[Nx-1][j][k];
	}
}

inline void LBM::Boundary(double delta, int bctype, Vector3d vw)
{
	// Vector3d v0 (0.01, 0., 0.); 
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		BounceBack(0 ,j,k);
		BounceBack(Nx,j,k);
	}

	for (int i=0; i<=Nx; ++i)
	for (int k=0; k<=Nz; ++k)
	{
		BounceBack(i,0 ,k);
		BounceBack(i,Ny,k);
	}

	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	{
		BounceBack(i,j,0 );
		BounceBack(i,j,Nz);
	}
}

inline void LBM::ApplyWall()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t c=0; c<Lwall.size(); ++c)
	{
		BounceBack(Lwall[c](0),Lwall[c](1),Lwall[c](2));
	}
}

inline void LBM::SetWall()
{
	// Vector3d v0 (0.01, 0., 0.); 
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		BounceBack(0 ,j,k);
		BounceBack(Nx,j,k);
		// F[0 ][j][k] = F[1   ][j][k];
		// F[Nx][j][k] = F[Nx-1][j][k];
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int k=0; k<=Nz; ++k)
	{
		BounceBack(i,0 ,k);
		BounceBack(i,Ny,k);
		// F[i][0 ][k] = F[i][1   ][k];
		// F[i][Ny][k] = F[i][Ny-1][k];
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	{
		BounceBack(i,j,0 );
		BounceBack(i,j,Nz);
		// F[i][j][0 ] = F[i][j][1   ];
		// F[i][j][Nz] = F[i][j][Nz-1];
	}
}

// Velocity interpolation based bounce back, make sure that 0<delta<1, q for the distribution that missed after streaming
inline void LBM::VIBB(int i, int j, int k, int q, double delta, Vector3d& vw, Vector3d& fh)
{
	int in = (i + (int) E[q][0]);
	int jn = (j + (int) E[q][1]);
	int kn = (k + (int) E[q][2]);

	// if (G[in][jn][kn]==-2.)
	// {
	// 	cout << "G[in][jn][kn]==-2." << endl;
	// 	cout << "i= " << i << " j= " << j << " k= " << k << endl;
	// }

	double fi, rhob;

	if (in<0 || in>Nx || jn<0 || jn>Ny || kn<0 || kn>Nz)
	{
		cout << "3213213213" << endl;
		fi = Ft[i][j][k](Op[q]);
		rhob = Rho[i][j][k];
	}

	else
	{
		double rhof 	= Rho[i][j][k];
		double rhoff 	= Rho[in][jn][kn];
		Vector3d vf 	= V[i][j][k];
		Vector3d vff 	= V[in][jn][kn];
		Vector3d vb 	= vw;
		Vector3d vi;
		if(delta<0.5) 	vi = 2*delta*vf + (1-2*delta)*vff;
		else			vi = ((1-delta)*vf + (2*delta-1)*vb)/delta;

		// Vector3d vii = (2*delta*vb + (1-delta)*vff)/(1+delta);
		// vi = (vi+2.*vii)/3.;

		double rhoi 	= 2*delta*rhof + (1-2*delta)*rhoff;
		rhob 	= delta*(rhof-rhoff) + rhof;

		VectorXd feqi(Q), feqf(Q), feqff(Q);
		(this->*CalFeq)(feqi, rhoi, vi);
		(this->*CalFeq)(feqf, rhof, vf);
		(this->*CalFeq)(feqff, rhoff, vff);

		// double fneqi	= 2*delta*(Ft[i][j][k](Op[q]) - feqf(Op[q]))/Rho[i][j][k]*rhoi + (1-2*delta)*(Ft[in][jn][kn](Op[q]) - feqff(Op[q]))/Rho[in][jn][kn]*rhoi;

		// double fneqiRatio	= 2*delta*(Ft[i][j][k](Op[q]) - feqf(Op[q]))/feqf(Op[q]) + (1-2*delta)*(Ft[in][jn][kn](Op[q]) - feqff(Op[q]))/feqff(Op[q]);

		double fneqi	= Ft[i][j][k](Op[q]) - feqf(Op[q])/*/Rho[i][j][k]*rhoi*/;
		
		// double fneqiRatio	= (Ft[i][j][k](Op[q]) - feqf(Op[q]))/feqf(Op[q]);

		fi = feqi(Op[q]) + fneqi;

		// if (i==42 && j==158)
		// {
		// 	cout << "vi= " << vi.transpose() << endl;
		// 	cout << "rhoi= " << rhoi << endl;
		// 	cout << "feqi(Op[q])= " << feqi(Op[q]) << endl;
		// 	cout << "fneqi= "<< fneqi << endl;
		// }
	}

	// double fi = feqi(Op[q])*(1.+fneqiRatio);

	F[i][j][k](q) = fi + 6*W[q]*rhob*E[q].dot(vw);

	if (F[i][j][k](q)<0.)
	{
		// F[i][j][k](q) = feqi(Op[q]) + (Ft[in][jn][kn](Op[q]) - feqff(Op[q])) + 6*W[q]*rhob*E[q].dot(vw);

		if (F[i][j][k](q)<0.)
		{
			cout << "i= " << i << " j= "<< j << " k= "<< k << endl;
			cout << "negative f in VIBB" << endl;
			cout << "fi= "<< fi << endl;
			// cout << "fneqff= "<< (Ft[in][jn][kn](Op[q]) - feqff(Op[q])) << endl;
			// cout << "feqi(Op[q])= "<< feqi(Op[q]) << endl;
			cout << "6*W[q]*rhob*E[q].dot(vw)= "<< 6*W[q]*rhob*E[q].dot(vw) << endl;

			abort();			
		}
		// cout << "negative f in VIBB" << endl;
		// cout << "fi= "<< fi << endl;
		// cout << "fneqi= "<< fneqi << endl;
		// cout << "feqi(Op[q])= "<< feqi(Op[q]) << endl;
		// cout << "6*W[q]*rhob*E[q].dot(vw)= "<< 6*W[q]*rhob*E[q].dot(vw) << endl;

		// abort();
	}

	// fh = -(F[i][j][k](q)+fi)*E[q];

	// fh = F[i][j][k](q)*(vw-E[q]) - fi*(E[q]+vw);

	// fh = F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);

	fh = -(F[i][j][k](q) + Ft[i][j][k](Op[q]))*E[q];

	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	// fh = F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);
	// cout << F[i][j][k](q) + Ft[i][j][k](Op[q]) << endl;
	// cout << "in vibb fh= " << fh.transpose() << endl;
}

inline void LBM::IBBYu(int i, int j, int k, int q, double delta, Vector3d& vw, Vector3d& fh)
{
	int in = (i + (int) E[q][0]);
	int jn = (j + (int) E[q][1]);
	int kn = (k + (int) E[q][2]);

	if (in<0 || in>Nx || jn<0 || jn>Ny || kn<0 || kn>Nz)
	{
		cout << "3213213213" << endl;
	}

	double fw = delta*Ft[i][j][k](Op[q]) + (1-delta)*Ft[in][jn][kn](Op[q]) + 6*W[q]*Rho[i][j][k]*E[q].dot(vw);
	F[i][j][k](q) = fw/(1+delta) + delta*F[in][jn][kn](q)/(1+delta);

	// Eq.42 (Peng 2016)
	// Implementation issues and benchmarking of lattice Boltzmann method for moving rigid particle simulations in a viscous flow
	// fh = F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);
	fh = -(F[i][j][k](q) + Ft[i][j][k](Op[q]))*E[q];
}

inline void LBM::PSM(int i, int j, int k, double g, Vector3d& vw, Vector3d& fh)
{
	double bn = g*(Tau-0.5)/(1.-g + Tau-0.5);
	// double bn = g;
	
	VectorXd feq(Q), feqv(Q);
	(this->*CalFeq)(feq, Rho[i][j][k], V[i][j][k]);
	(this->*CalFeq)(feqv, Rho[i][j][k], vw);

	for (int q=0; q<Q; ++q)
	{
		double omegas = F[i][j][k](Op[q])- F[i][j][k](q) + feqv(q) - feq(Op[q]);
		Ft[i][j][k](q) = F[i][j][k](q) + (1.-bn)*Omega0*(feq(q) - F[i][j][k](q)) + bn*omegas;
		fh += omegas*E[q];
	}
	fh *= -bn;
}

inline void LBM::VAM(int i, int j, int k, double g, double rhos, Vector3d& vs, Vector3d& fh)
{
	// volume average velocity
	Vector3d vm = (g*rhos*vs + (1-g)*Rho[i][j][k]*V[i][j][k])/(g*rhos + (1-g)*Rho[i][j][k]);
	// Vector3d vm = g*vs + (1-g)*V[i][j][k];
	// double re = (V[i][j][k]-vm).norm()/Nu
	// momentum transfer
	fh = (V[i][j][k]-vm)*Rho[i][j][k]*(1.-g);
	// apply though body force
	BodyForceLocalOpenMP(i, j, k, -fh);
	// #pragma omp critical
	// {
	//     BodyForceLocal(i, j, k, -fh);
	//     V[i][j][k] = vm;
	// }
}

// Linear interpolate velocity at any poit x
inline Vector3d LBM::InterpolateV(const Vector3d& x)
{
    Vector3d vx (0., 0., 0.);
    // Linear interpolation
    // Find node
    Vector3i min0, max0;
    min0(0) = floor(x(0));
    max0(0) = min0(0)+1;
    min0(1) = floor(x(1));
    max0(1) = min0(1)+1;
    min0(2) = floor(x(2));
    max0(2) = min0(2)+1;

    // for(int d=0; d<D; ++d)
    // {
    // 	if (Periodic[d])
    // 	{
    // 		min0(d) = (min0(d)+DomSize[d]+1)%(DomSize[d]+1);
    // 		max0(d) = (max0(d)+DomSize[d]+1)%(DomSize[d]+1);
    // 	}
    // }

    vector <Vector3d> ver;
    ver.resize(8);

    ver[0] << min0(0), min0(1), min0(2);
    ver[1] << max0(0), min0(1), min0(2);
    ver[2] << max0(0), max0(1), min0(2);
    ver[3] << min0(0), max0(1), min0(2);

    ver[4] << max0(0), min0(1), max0(2);
    ver[5] << min0(0), min0(1), max0(2);
    ver[6] << min0(0), max0(1), max0(2);
    ver[7] << max0(0), max0(1), max0(2);

    for (size_t l=0; l<pow(2,D); ++l)
    {
        int i = (int) ver[l](0);
        int j = (int) ver[l](1);
        int k = (int) ver[l](2);
        i = (i+Nx+1)%(Nx+1);
        j = (j+Ny+1)%(Ny+1);
        k = (k+Nz+1)%(Nz+1);
        Vector3d s = x-ver[7-l];
        double vol = abs(s(0)*s(1)*s(2));
        vx += vol*V[i][j][k];
    }
    // if (vx.norm()>0.05)
    // {
    // 	cout << "x: " << x.transpose() << endl;
    // 	cout << "vx: " << vx.transpose() << endl;
    // 	for (size_t l=0; l<pow(2,D); ++l)
    // 	{
    // 		int i = (int) ver[l](0);
    //     	int j = (int) ver[l](1);
    //     	int k = (int) ver[l](2);
    //     	Vector3d s = x-ver[7-l];
    //     	cout << "s: " << s.transpose() << endl;
    //     	double vol = abs(s(0)*s(1)*s(2));
    //     	cout << "vol: " << vol << endl;
    // 		cout << "ver[l]: " << ver[l].transpose() << endl;	
    // 		cout << "V[i][j][k]: " << V[i][j][k].transpose() << endl;	
    // 	}
    // 	abort();
    // }
    return vx;
}

// inline double LBM::Bn0(double g)
// {
// 	return g;
// }

// inline double LBM::Bn1(double g)
// {
// 	return g*(Tau-0.5)/(1.-g + Tau-0.5);
// }

/*inline void LBM::Refill(int i, int j, int k, Vector3d vw, Vector3d fh)
{

}*/

// inline void LBM::WriteFileH5(int n)
// {
// 	stringstream	out;					//convert int to string for file name.
// 	out << setw(6) << setfill('0') << n;			
// 	string file_name_h5 = "LBM"+out.str()+".h5";

// 	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
// 	int numLat = (Nx+1)*(Ny+1)*(Nz+1);

// 	hsize_t	dims_scalar[3] = {Nz+1, Ny+1, Nx+1};			//create data space.

// 	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);

// 	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);

// 	double* rho_h5 	= new double[numLat];
// 	double* g_h5 	= new double[numLat];
// 	double* vel_h5 	= new double[numLat];
// 	double* vx_h5 	= new double[numLat];
// 	double* vy_h5 	= new double[numLat];
// 	double* vz_h5 	= new double[numLat];

// 	int len = 0;
// 	// #pragma omp parallel for schedule(static) num_threads(Nproc)
// 	for (int k=0; k<=Nz; k++)
// 	for (int j=0; j<=Ny; j++)
// 	for (int i=0; i<=Nx; i++)
// 	{
//         rho_h5[len] = Rho[i][j][k]; 
//         // g_h5[len] = G[i][j][k][1] + G[i][j][k][3];
//         g_h5[len] = G[i][j][k][0];

// 		vel_h5[len] = V[i][j][k].norm();

// 		vx_h5[len] = V[i][j][k](0);
// 		vy_h5[len] = V[i][j][k](1);
// 		vz_h5[len] = V[i][j][k](2);
//         len++;
// 	}

// 	DataSet	*dataset_rho = new DataSet(file.createDataSet("Density", PredType::NATIVE_DOUBLE, *space_scalar));
// 	DataSet	*dataset_g	 = new DataSet(file.createDataSet("Gamma", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_vx  = new DataSet(file.createDataSet("Vx", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_vy  = new DataSet(file.createDataSet("Vy", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_vz  = new DataSet(file.createDataSet("Vz", PredType::NATIVE_DOUBLE, *space_scalar));

// 	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
// 	dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
// 	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
// 	dataset_vx->write(vx_h5, PredType::NATIVE_DOUBLE);
// 	dataset_vy->write(vy_h5, PredType::NATIVE_DOUBLE);
// 	dataset_vz->write(vz_h5, PredType::NATIVE_DOUBLE);

// 	delete space_scalar;
// 	delete dataset_rho;
// 	delete dataset_g;
// 	delete dataset_vel;
// 	delete dataset_vx;
// 	delete dataset_vy;
// 	delete dataset_vz;
// 	delete rho_h5;
// 	delete g_h5;
// 	delete vel_h5;
// 	delete vx_h5;
// 	delete vy_h5;
// 	delete vz_h5;

// 	file.close();

// 	string file_name_xmf = "LBM_"+out.str()+".xmf";

//     std::ofstream oss;
//     oss.open(file_name_xmf);
//     oss << "<?xml version=\"1.0\" ?>\n";
//     oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
//     oss << "<Xdmf Version=\"2.0\">\n";
//     oss << " <Domain>\n";
//     oss << "   <Grid Name=\"LBM\" GridType=\"Uniform\">\n";
//     oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\"/>\n";
//     oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
//     oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
//     oss << "       </DataItem>\n";
//     oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0  << "\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Geometry>\n";
//     oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Density \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Velocity \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Gamma \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "   </Grid>\n";
//     oss << " </Domain>\n";
//     oss << "</Xdmf>\n";
//     oss.close();
// }

inline void LBM::WriteFileH5(int n, int scale)
{
	stringstream	out;					//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "LBM"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.

	hsize_t nx = (Nx+1)/scale;
	hsize_t ny = (Ny+1)/scale;
	hsize_t nz = (Nz+1)/scale;

	int numLat = nx*ny*nz;

	hsize_t	dims_scalar[3] = {nx, ny, nz};			//create data space.
	hsize_t	dims_vector[4] = {nx, ny, nz, 3};		//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* rho_h5 	= new double[numLat];
	double* g_h5 	= new double[numLat];
	double* vel_h5 	= new double[3*numLat];

	int len = 0;
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t k=0; k<nz; k++)
	for (size_t j=0; j<ny; j++)
	for (size_t i=0; i<nx; i++)
	{		
		rho_h5[len] = 0.;
		g_h5[len] = 0.;
		vel_h5[3*len  ] = 0.;
		vel_h5[3*len+1] = 0.;
		vel_h5[3*len+2] = 0.;
		int cout = 0;
		for (int kk=0; kk<scale; kk++)
		for (int jj=0; jj<scale; jj++)
		for (int ii=0; ii<scale; ii++)
		{
			if (i+ii<=(size_t) Nx && j+jj<=(size_t) Ny && k+kk<=(size_t) Nz)
			{
				rho_h5[len] += Rho[i+ii][j+jj][k+kk];
				g_h5[len] += G[i+ii][j+jj][k+kk][0];
				vel_h5[3*len  ] += V[i][j][k](0);
				vel_h5[3*len+1] += V[i][j][k](1);
				vel_h5[3*len+2] += V[i][j][k](2);
				cout++;
			}
		}
		rho_h5[len] /= (double) cout;
		g_h5[len] /= (double) cout;
		vel_h5[3*len  ] /= (double) cout;
		vel_h5[3*len+1] /= (double) cout;
		vel_h5[3*len+2] /= (double) cout;
        len++;
	}

	DataSet	*dataset_rho = new DataSet(file.createDataSet("Density", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_g	 = new DataSet(file.createDataSet("Gamma", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete dataset_rho;
	delete dataset_g;
	delete dataset_vel;

	delete rho_h5;
	delete g_h5;
	delete vel_h5;

	file.close();

	string file_name_xmf = "LBM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"LBM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << nz << " " << ny << " " << nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Density \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Gamma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}