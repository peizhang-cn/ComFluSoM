#include "../HEADER.h"
#include <LBM.h>

class MLBM
{
public:
	MLBM();
	~MLBM();
	MLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, double dc);
	void Init(double rho0, double c0, Vector3d initV);
	void Solve(int tt, int ts, int scale);
	void WriteFileH5(int n, int scale);

	LBM*				DomF;				// Domain of fluid
	LBM*				DomS;				// Domain of scalar	
};

inline MLBM::MLBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, double dc)
{
	DomF = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, nu);
	DomS = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, dc);
}

inline void MLBM::Init(double rho0, double c0, Vector3d initV)
{
	DomF->Init(rho0, initV);
	DomS->InitCDE(c0, initV);
	cout << "finish InitCDE" << endl;
	DomS->V = DomF->V;
}

inline void MLBM::Solve(int tt, int ts, int scale)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts==0)
		{
			cout << "Time= " << t << endl;
			WriteFileH5(t,scale);			
		}
		DomF->CollideSRT();
		DomS->CollideSRTCDE();
		DomF->ApplyWall();
		DomS->ApplyWall();
		DomF->Stream();
		DomS->Stream();
		DomF->ApplyNoGradBC();
		DomF->CalRhoV();	
		DomS->CalRho();
		DomF->ApplyVelBC();
	}
}

inline void MLBM::WriteFileH5(int n, int scale)
{
	stringstream	out;					//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "MLBM"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.

	hsize_t nx = (DomF->Nx+1)/scale;
	hsize_t ny = (DomF->Ny+1)/scale;
	hsize_t nz = (DomF->Nz+1)/scale;

	int numLat = nx*ny*nz;

	hsize_t	dims_scalar[3] = {nx, ny, nz};			//create data space.
	hsize_t	dims_vector[4] = {nx, ny, nz, 3};		//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* rho_h5 	= new double[numLat];
	double* c_h5 	= new double[numLat];
	double* g_h5 	= new double[numLat];
	double* vel_h5 	= new double[3*numLat];

	int len = 0;
	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t k=0; k<nz; k++)
	for (size_t j=0; j<ny; j++)
	for (size_t i=0; i<nx; i++)
	{		
		rho_h5[len] = 0.;
		c_h5[len] = 0.;
		g_h5[len] = 0.;
		vel_h5[3*len  ] = 0.;
		vel_h5[3*len+1] = 0.;
		vel_h5[3*len+2] = 0.;
		int cout = 0;
		for (int kk=0; kk<scale; kk++)
		for (int jj=0; jj<scale; jj++)
		for (int ii=0; ii<scale; ii++)
		{
			if (i+ii<=(size_t) DomF->Nx && j+jj<=(size_t) DomF->Ny && k+kk<=(size_t) DomF->Nz)
			{
				rho_h5[len] += DomF->Rho[i+ii][j+jj][k+kk];
				c_h5[len] += DomS->Rho[i+ii][j+jj][k+kk];
				g_h5[len] += DomF->G[i+ii][j+jj][k+kk][0];
				vel_h5[3*len  ] += DomF->V[i][j][k](0);
				vel_h5[3*len+1] += DomF->V[i][j][k](1);
				vel_h5[3*len+2] += DomF->V[i][j][k](2);
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
	DataSet	*dataset_c 	 = new DataSet(file.createDataSet("Concentration", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_g	 = new DataSet(file.createDataSet("Gamma", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
	dataset_c->write(c_h5, PredType::NATIVE_DOUBLE);
	dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete dataset_rho;
	delete dataset_c;
	delete dataset_g;
	delete dataset_vel;

	delete rho_h5;
	delete c_h5;
	delete g_h5;
	delete vel_h5;

	file.close();

	string file_name_xmf = "MLBM_"+out.str()+".xmf";

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
    oss << "     <Attribute Name=\"Concentration\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Concentration \n";
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