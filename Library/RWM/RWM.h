#include "../HEADER.h"
#include <RWM_PARTICLE.h>

class RWM
{
public:

//Public Variables =====================================================================================================================================

	bool 							Periodic[3];
	int             			   	Nproc;                                                   	// Number of processors which used
    vector <RWM_PARTICLE*>         	Lp;                                                      	// List of RWM particles
    double***                      	C;                                                     		// Scalar concentration
    Vector3d***				   		V_ptr;														// Pointer to the velocity field

//Public Functions ======================================================================================================================================

	RWM(): gen(std::random_device()()) {};
	RWM(int nx, int ny, int nz, double dc, double dt);
	double GetNormalD();
	double GetUniformD0();
	double GetUniformD1();
	Vector3d GetRandomPointinUnitSphere();
	void Init();
	void Move();
	void CalV(const Vector3d& x, Vector3d*** v_ptr, Vector3d& vx);
	Vector3d InterpolateV(const Vector3d& x);
	void CalC();
	double CalP(double k);
	void RemoveParticles();
	void ContactPointWithSphere(RWM_PARTICLE* p0, Vector3d xs/*centre of sphere*/, Vector3d xsb, double r/*radius*/, Vector3d& xc, Vector3d& n);
	void Reflection(RWM_PARTICLE* p0, Vector3d& n/*normal vector*/, Vector3d& xc/*contact point*/);
	void Reaction(RWM_PARTICLE* p0, double Possibility, Vector3d& n/*normal vector*/, Vector3d& xc/*contact point*/);
	void AddParticle(int tag, Vector3d& x, double m);
	void WriteFileH5(int n);
	
// private:

//Private Variables =====================================================================================================================================

    int                            	Nx;                                                      	// Domain size
    int                            	Ny;
    int                            	Nz;
    vector <int>					DomSize;													// Domain size for loop
    int                            	D;                                                       	// Dimension
    double							Dc;															// Diffusion coefficient
    double 							Dt;															// Time step
	double 							K1th;														// First order reaction coefficient
    double 							P;															// Possibility of reaction	
	mt19937 					   	gen;														// Random number generator
	vector <int>					Lrm;														// List for removed particles
};

// Standard normal distribution
inline double RWM::GetNormalD()
{
	normal_distribution<double>			dis(0.,1.);
	return dis(gen);
}

// Uniform distribution from 0 to 1
inline double RWM::GetUniformD0()
{
	uniform_real_distribution<double>	dis(0.,1.);
	return dis(gen);
}

// Uniform distribution from -0.5 to 0.5
inline double RWM::GetUniformD1()
{
	uniform_real_distribution<double>	dis(-0.5,0.5);
	return dis(gen);
}
// get a random point in a unit sphere (circle in 2D)
inline Vector3d RWM::GetRandomPointinUnitSphere()
{
	double u = GetUniformD0();
	Vector3d x = Vector3d::Zero();
	for (size_t d=0; d<D; ++d)
	{
		x(d) = GetNormalD();
	}
	x.normalized();
	x *= pow(u, 1/D);
	return x;
}

inline RWM::RWM(int nx, int ny, int nz, double dc, double dt) : gen(std::random_device()())
{
	Nx = nx;
	Ny = ny;
	Nz = nz;
	D=3;
	if (Nz==0)
	{
		D=2;
		if (Ny==0)	D=1;
	}

	Nproc 	= 1;
	Dc 		= dc;
	Dt 		= dt;

	Lp.resize(0);

	Periodic[0] = true;
	Periodic[1] = true;
	Periodic[2] = true;
}

inline void RWM::Init()
{
	DomSize.push_back(Nx);
	DomSize.push_back(Ny);
	DomSize.push_back(Nz);

	Lrm.resize(0);

	// Init concentration field
	C = new double** [Nx+1];
	for (int i=0; i<=Nx; ++i)
	{
		C[i] = new double* [Ny+1];
		for (int j=0; j<=Ny; ++j)
		{
			C[i][j]	= new double [Nz+1];
			for (int k=0; k<=Nz; ++k)
			{
				C[i][j][k] = 0.;
			}
		}
	}
}

inline void RWM::Move()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		// Store the position before move
		Lp[p]->Xb = Lp[p]->X;

		// Vector3d vp (0.,0.,0.);
		// cout << "start CalV" << endl;
		// CalV(Lp[p]->X, V_ptr, vp);

		// convection velocity
		Vector3d vc = InterpolateV(Lp[p]->X);

		for(int d=0; d<D; ++d)
		{
			Lp[p]->X(d) += vc(d)*Dt + GetNormalD()*sqrt(2*Dc*Dt);
			// if (Lp[p]->X(d)<0 || Lp[p]->X(d)>DomSize[d])	Lrm.push_back(p);
		}
		Lp[p]->Vc.setZero();
	}
}

inline void RWM::ContactPointWithSphere(RWM_PARTICLE* p0, Vector3d xs/*centre of sphere*/, Vector3d xsb, double r/*radius*/, Vector3d& xc, Vector3d& n)
{
	Vector3d aa = p0->X - p0->Xb - (xs-xsb);
	Vector3d ac = p0->Xb - xsb;
	double A = aa.squaredNorm();
	double B = 2*ac.dot(aa);
	double C = ac.squaredNorm()-r*r;
	double deltaRoot = sqrt(B*B-4*A*C);
	double k0 = (-B+deltaRoot)/(2*A);
	double k1 = (-B-deltaRoot)/(2*A);
	double k = k0;
	if (k1>=0.&&k1<=1.)		k=k1;
	if (k<0.|| k>1.)
	{
		double kmin = min(k0,k1);
		double kmax = max(k0,k1);
		if (abs(kmin)<abs(kmax-1))	k = 0.;
		else 						k = 1.;
		// cout << "wrong k: " << k << endl;
		// cout << "k1: " << k0 << endl;
		// cout << "k2: " << k1 << endl;
		// cout << "x: " << p0->X.transpose() << endl;
		// cout << "xb: " << p0->Xb.transpose() << endl;
		// cout << "xs: " << xs.transpose() << endl;
		// cout << "xsb: " << xsb.transpose() << endl;
		// cout << (p0->Xb-xsb).norm() << endl;
		// cout << (p0->X-xs).norm() << endl;
		// abort();
	}
	xc = p0->Xb + k*(p0->X - p0->Xb);
	n = (xc-xsb-k*(xs-xsb)).normalized();
}

// Make true that the RW particle contact with walls before use this function
inline void RWM::Reflection(RWM_PARTICLE* p0, Vector3d& n/*normal vector*/, Vector3d& xc/*contact point*/)
{
	// http://mathworld.wolfram.com/Reflection.html
	// Refected vector
	Vector3d v = -p0->Xb+xc + 2.*n.dot(p0->Xb-xc)*n;
	// Reflected position
	p0->X = (p0->X-xc).norm()*v.normalized()+xc;
}

inline void RWM::Reaction(RWM_PARTICLE* p0, double Possibility, Vector3d& n/*normal vector*/, Vector3d& xc/*contact point*/)
{
	if (GetUniformD0()>Possibility)		Lrm.push_back(p0->ID);
	else								Reflection(p0, n, xc);
}

// An improved scheme for a Robin boundary condition in discrete-time random walk algorithms
inline double RWM::CalP(double k /*1th order recation coefficient*/)
{
	double p1 = k*sqrt(M_PI*Dt/Dc);
	return (p1/(1.+0.5*p1));
}

inline void RWM::RemoveParticles()
{
	// Sort Lrm from large to small.
	sort(Lrm.begin(), Lrm.end(), greater<int>());
	Lrm.erase( unique( Lrm.begin(), Lrm.end() ), Lrm.end() );
	// Switch removed particle with last particle then remove it.
	for (size_t p=0; p<Lrm.size(); ++p)
	{
		Lp[Lrm[p]] = Lp.back();
		// delete Lp.back();
		Lp.pop_back();	
	}
	// Reset Lrm to zero.
	Lrm.resize(0);
}

inline void RWM::CalV(const Vector3d& x, Vector3d*** v_ptr, Vector3d& vx)
{
	for(int d=0; d<D; ++d)
	{
		if (x(d)<0 || x(d)>DomSize[d])
		{
			cout << "CalV out of domain" << endl;
		}
	}
	vx << 0., 0., 0.;
	
	// Linear interpolation
	// Find node
	int xmin = int(x(0));
	int xmax = xmin+1;
	int ymin = int(x(1));
	int ymax = ymin+1;
	int zmin = int(x(2));
	int zmax = zmin+1;

	cout << "xmin " << xmin << endl;
	cout << "xmax " << xmax << endl;
	cout << "ymin " << ymin << endl;
	cout << "ymax " << ymax << endl;
	cout << "zmin " << zmin << endl;
	cout << "zmax " << zmax << endl;

	if (Periodic[0])
	{
		xmin = (xmin+Nx+1)%(Nx+1);
		xmax = (xmax+Nx+1)%(Nx+1);
	}
	if (Periodic[1])
	{
		ymin = (ymin+Ny+1)%(Ny+1);
		ymax = (ymax+Ny+1)%(Ny+1);
	}
	// if (Periodic[2])
	// {
	// 	zmin = (zmin+Nz+1)%(Nz+1);
	// 	zmax = (zmax+Nz+1)%(Nz+1);
	// }

	vector <Vector3d> ver;

	ver.resize(8);

	ver[0] << xmin, ymin, zmin;
	ver[1] << xmax, ymin, zmin;
	ver[2] << xmax, ymax, zmin;
	ver[3] << xmin, ymax, zmin;

	ver[4] << xmax, ymin, zmax;
	ver[5] << xmin, ymin, zmax;
	ver[6] << xmin, ymax, zmax;
	ver[7] << xmax, ymax, zmax;

	cout << "zmax " << zmax << endl;
	cout << "ver[7] " << ver[5].transpose() << endl;

	cout << "start 1111" << endl;
	cout << "x " << x.transpose() << endl;
	for (size_t l=0; l<pow(2,D); ++l)
	{
		int i = (int) ver[l](0);
		int j = (int) ver[l](1);
		int k = (int) ver[l](2);
		cout << "l " << l << endl;
		Vector3d s = x-ver[7-l];
		double vol = abs(s(0)*s(1)*s(2));
		cout << "vol " << vol << endl;
		cout << i << " " << j << " " << k << endl;
		cout << "s " << s.transpose() << endl;
		cout << "ver[7-l] " << ver[7-l].transpose() << endl;
		cout << "ver[7] " << ver[7].transpose() << endl;
		// abort();
		cout << v_ptr[i][j][k] << endl;
		vx += vol*(v_ptr[i][j][k]);
	}
	abort();
}

// Linear interpolate velocity at any poit x
inline Vector3d RWM::InterpolateV(const Vector3d& x)
{
    Vector3d vx (0., 0., 0.);
    // Linear interpolation
    // Find node
    Vector3i min0, max0;
    min0(0) = int(x(0));
    max0(0) = min0(0)+1;
    min0(1) = int(x(1));
    max0(1) = min0(1)+1;
    min0(2) = int(x(2));
    max0(2) = min0(2)+1;

    for(int d=0; d<D; ++d)
    {
    	if (Periodic[d])
    	{
    		min0(d) = (min0(d)+DomSize[d]+1)%(DomSize[d]+1);
    		max0(d) = (max0(d)+DomSize[d]+1)%(DomSize[d]+1);
    	}
    }

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
        Vector3d s = x-ver[7-l];
        double vol = abs(s(0)*s(1)*s(2));
        vx += vol*V_ptr[i][j][k];
    }
    return vx;
}

inline void RWM::CalC()
{
	// Reset concentration to zero
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		C[i][j][k] = 0.;
	}

	// Sum particle mass to concentration
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		int i = round(Lp[p]->X(0));
		int j = round(Lp[p]->X(1));
		int k = round(Lp[p]->X(2));
		i = (i+Nx+1)%(Nx+1);
		j = (j+Ny+1)%(Ny+1);
		k = (k+Nz+1)%(Nz+1);
		#pragma omp atomic
		C[i][j][k] += Lp[p]->M;
	}

	// For boundary nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int j=0; j<=Ny; ++j)
	{
		C[i][j][0 ] /= 2.;
		C[i][j][Nz] /= 2.;
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<=Nx; ++i)
	for (int k=0; k<=Nz; ++k)
	{
		C[i][0 ][k] /= 2.;
		C[i][Ny][k] /= 2.;
	}
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int j=0; j<=Ny; ++j)
	for (int k=0; k<=Nz; ++k)
	{
		C[0 ][j][k] /= 2.;
		C[Nx][j][k] /= 2.;
	}
}

inline void RWM::AddParticle(int tag, Vector3d& x, double m)
{
    Lp.push_back(new RWM_PARTICLE(tag,x,m));
    Lp[Lp.size()-1]->ID = Lp.size()-1;
}

inline void RWM::WriteFileH5(int n)
{
	stringstream	out;					//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "RWM"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	int numLat = (Nx+1)*(Ny+1)*(Nz+1);

	hsize_t	dims_scalar[3] = {Nz+1, Ny+1, Nx+1};			//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);

	double* c_h5 	= new double[numLat];

	int len = 0;
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int k=0; k<=Nz; k++)
	for (int j=0; j<=Ny; j++)
	for (int i=0; i<=Nx; i++)
	{
        c_h5[len] = C[i][j][k];
        len++;
	}

	DataSet	*dataset_c = new DataSet(file.createDataSet("Concentration", PredType::NATIVE_DOUBLE, *space_scalar));

	dataset_c->write(c_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete dataset_c;

	delete c_h5;

	file.close();

	string file_name_xmf = "RWM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"RWM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Concentration\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz+1 << " " << Ny+1 << " " << Nx+1 << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Concentration \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}