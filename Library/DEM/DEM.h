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
	DEM(int nx, int ny, int nz);
	~DEM();
	void Init();
	void AddSphere(int tag, double r, Vector3d& x, double rho);
	void AddDisk2D(int tag, double r, Vector3d& x, double rho);
	void Move(bool first, double dt);
	void ZeroForceTorque();
	void SetG(Vector3d& g);
	void RecordX();																			// Record position at Xb for check refilling LBM nodes
	void Contact2P(DEM_PARTICLE* p0, DEM_PARTICLE* p1);
	void FindContact();
	void Contact();
	void Solve(int tt, int ts, double dt);
	void WriteFileH5(int n);

	int 							Nproc;
    int 							Nx;														// Mesh size for contact detection
    int 							Ny;
    int 							Nz;

    int***		 					Flag;													// Flag of lattice type

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < vector< int > >        Lc;                                                 	// List of potential contacted paricles' ID
};

inline DEM::DEM(int nx, int ny, int nz)
{
	Nx = nx;
	Ny = ny;
	Nz = nz;

	Nproc = 1;
	// Np	= 0;
}

inline void DEM::Init()
{
    cout << "================ Start init. ================" << endl;

	Flag = new int** [Nx+1];

	for (int i=0; i<=Nx; ++i)
	{
		Flag[i] = new int* [Ny+1];

		for (int j=0; j<=Ny; ++j)
		{
			Flag[i][j]	= new int [Nz+1];

			for (int k=0; k<=Nz; ++k)
			{
				Flag[i][j][k] = -1;
			}
		}
	}
	cout << "================ Finish init. ================" << endl;
}

inline void DEM::AddSphere(int tag, double r, Vector3d& x, double rho)
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		if (tag==Lp[i]->Tag)
		{
			cout << "\033[1;31mError: Cannot add sphere. The tag is existed.\033[0m\n";		
			exit(0);
		}
	}
	Lp.push_back(new SPHERE(tag, r, x, rho));
	Lp[Lp.size()-1]->ID = Lp.size()-1;
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

inline void DEM::Move(bool firstStep, double dt)
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->VelocityVerlet(firstStep, dt);
	}
}

inline void DEM::ZeroForceTorque()
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->ZeroForceTorque();
	}
}

inline void DEM::SetG(Vector3d& g)
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->SetG(g);
	}
}

inline void DEM::RecordX()
{
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->Xb = Lp[i]->X;
	}
}

// Contact force model, only normal force is considered now.
inline void DEM::Contact2P(DEM_PARTICLE* p0, DEM_PARTICLE* p1)
{
	// Normal direction (p0 pinnts to p1)
	Vector3d n = p1->X-p0->X;
	// Overlapping distance
	double delta = p0->R+p1->R-n.norm();

	if (delta>0.)
	{
		n.normalize();
		Vector3d fc = n*pow(delta, 1.5);
		p0->Fc -= p0->Kn*fc;
		p1->Fc += p1->Kn*fc;
	}
}

inline void DEM::FindContact()
{
	for (int p=0; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];

		if (p0->MinX<0. || p0->MaxX>Nx || p0->MinY<0. || p0->MaxY>Ny || p0->MinZ<0. || p0->MaxZ>Nz)
		{
			cout << "\033[1;31mError: particle out of box.\033[0m\n";		
			exit(0);
		}

		for (int i=p0->MinX; i<=p0->MaxX; ++i)
		for (int j=p0->MinY; j<=p0->MaxY; ++j)
		for (int k=p0->MinZ; k<=p0->MaxZ; ++k)
		{
			Flag[i][j][k] = -1;
		}
	}

	for (int p=0; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];

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
				if (Flag[i][j][k]==-1)
				{
					Flag[i][j][k] = p0->ID;
				}
				else
				{
					int q = Flag[i][j][k];
					Lc.push_back({min(p,q), max(p,q)});

					if (p==q)
					{
						cout << "i= " << i << endl;
						cout << "j= " << j << endl;
						cout << "k= " << k << endl;

						cout << "p= " << p << endl;
						cout << "q= " << q << endl;
						abort();
					}
				}
			}
		}
	}

    sort( Lc.begin(), Lc.end() );
    Lc.erase( unique( Lc.begin(), Lc.end() ), Lc.end() );
}

inline void DEM::Contact()
{
    if (Lc.size()>0.)
    {
		for (size_t l=0; l<Lc.size(); ++l)
		{
			Contact2P(Lp[Lc[l][0]], Lp[Lc[l][1]]);
		}
    }
	Lc.clear();
}

inline void DEM::Solve(int tt, int ts, double dt)
{
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			WriteFileH5(t);
		}
		FindContact();
		Contact();
		if (t==0)	Move(true, dt);
		else 		Move(false, dt);
		ZeroForceTorque();
	}
}

inline void DEM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "DEM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};					//create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};				//create data space.

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

	double* r_h5 	= new double[  Lp.size()];
	double* tag_h5 	= new double[  Lp.size()];
	double* pos_h5 	= new double[3*Lp.size()];
	double* vel_h5 	= new double[3*Lp.size()];
	double* agv_h5	= new double[3*Lp.size()];

	for (size_t i=0; i<Lp.size(); ++i)
	{
		Vector3d agv = Lp[i]->Q._transformVector(Lp[i]->W);
        r_h5  [  i  ] 	= Lp[i]->R;
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
	}

	DataSet	*dataset_r 		= new DataSet(file.createDataSet("Radius", PredType::NATIVE_DOUBLE, *space_scalar));
	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_agv	= new DataSet(file.createDataSet("AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector));

	dataset_r->write(r_h5, PredType::NATIVE_DOUBLE);
	dataset_tag->write(tag_h5, PredType::NATIVE_DOUBLE);
	dataset_pos->write(pos_h5, PredType::NATIVE_DOUBLE);
	dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);
	dataset_agv->write(agv_h5, PredType::NATIVE_DOUBLE);

	delete space_scalar;
	delete space_vector;
	delete dataset_r;
	delete dataset_tag;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_agv;

	delete r_h5;
	delete tag_h5;
	delete pos_h5;
	delete vel_h5;
	delete agv_h5;

	file.close();

	string file_name_xmf = "DEM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"DEM\" GridType=\"Uniform\">\n";
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