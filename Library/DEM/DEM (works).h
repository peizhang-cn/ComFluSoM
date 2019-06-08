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
	void Solve(int tt, int ts, double dt);
	void WriteFileH5(int n);
	void WriteFileParticleInfo(int n);

	bool 							Periodic[3];

	int 							Nproc;
    int 							Nx;														// Mesh size for contact detection
    int 							Ny;
    int 							Nz;

    int 							DomSize[3];

    double 							Dt;														// Time step

    int***		 					Flag;													// Flag of lattice type

	vector < DEM_PARTICLE* >		Lp;														// List of particles
	vector < vector< size_t > >     Lc;                                                 	// List of potential contacted paricles' ID
	
	unordered_map<size_t, bool> 	CMap;													// Contact Map
	unordered_map<size_t, Vector3d> FMap;													// Friction Map
};

inline DEM::DEM(int nx, int ny, int nz)
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
	// Mark flag and add particles for the surounding box to handle walls and periodic BC
	// For x axis
	if (!Periodic[0])
	{
		for (int j=0; j<=Ny; ++j)
		for (int k=0; k<=Nz; ++k)
		{
			Flag[0 ][j][k] = 0;
			Flag[Nx][j][k] = 1;
		}
	}
	Vector3d x0 (0., 0., 0.);
	Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
	Lp[0]->ID = 0;
	Lp[0]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
	Lp[1]->ID = 1;
	Lp[1]->Type = -1;
	// For y axis
	if (!Periodic[1])
	{
		for (int i=0; i<=Nx; ++i)
		for (int k=0; k<=Nz; ++k)
		{
			Flag[i][0 ][k] = 2;
			Flag[i][Ny][k] = 3;
		}		
	}
	Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
	Lp[2]->ID = 2;
	Lp[2]->Type = -1;
	Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
	Lp[3]->ID = 3;
	Lp[3]->Type = -1;
	// For z axis
	if (!Periodic[2])
	{
		for (int i=0; i<=Nx; ++i)
		for (int j=0; j<=Ny; ++j)
		{
			Flag[i][j][0 ] = 4;
			Flag[i][j][Nz] = 5;
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
	// for (size_t i=0; i<Lp.size(); ++i)
	// {
	// 	if (tag==Lp[i]->Tag)
	// 	{
	// 		cout << "\033[1;31mError: Cannot add sphere. The tag is existed.\033[0m\n";		
	// 		exit(0);
	// 	}
	// }
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
		else if (p0->X(1)<1.)	p0->X(1) = p0->X(1)+Ny+1;
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
	for (size_t i=0; i<Lp.size(); ++i)
	{
		Lp[i]->Xb = Lp[i]->X;
	}
}

// Contact force model for spheres
inline void DEM::Contact2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, Vector3d& xi, bool& contacted)
{
	contacted = false;
	if (pi->Type==-1)
	{
		pi->X = pj->X;
		int dirc = pi->ID%2;
		int axis = (pi->ID-dirc)/2;
		if (!Periodic[axis])	pi->X(axis) = dirc*DomSize[axis];
		else 					pi->X(axis) += pj->R+1.;
	}
	// Normal direction (pj pinnts to pi)
	Vector3d n = pi->X-pj->X;
	// Overlapping distance
	double delta = pi->R+pj->R-n.norm();

	if (delta<0. && pi->Type!=-1)
	{
		Vector3d n2 = pi->X-pj->Xmir;
		double delta2 = pi->R+pj->R-n2.norm();
		Vector3d n3 = pi->Xmir-pj->X;
		double delta3 = pi->R+pj->R-n3.norm();
		if (delta2>0)
		{
			n = n2;
			delta = delta2;
		}
		if (delta3>0)
		{
			n = n3;
			delta = delta3;
		}
	}

	if (delta>0.)
	{
		if (delta>0.1)
		{
			cout << "delta= " << delta << endl;
			cout << pi->ID << endl;
			cout << pj->ID << endl;
			cout << pi->X.transpose() << endl;
			cout << pj->X.transpose() << endl;
			cout << pi->Xmir.transpose() << endl;
			cout << pj->Xmir.transpose() << endl;
			abort();			
		}
		// report contact for updating friction map 
		contacted = true;
		n.normalize();
		double kn 	= EffectiveValue(pi->Kn, pj->Kn);
		double kt 	= EffectiveValue(pi->Kt, pj->Kt);
		double gn 	= EffectiveValue(pi->Gn, pj->Gn);
		double gt 	= EffectiveValue(pi->Gt, pj->Gt);
		double me 	= EffectiveValue(pi->M, pj->M);
		// For collision with wall
		if (pi->Type==-1) 	me = pj->M;
 		// Relative velocity in normal direction
		Vector3d vn = (pj->V-pi->V).dot(n)*n;
		// Normal contact force
		Vector3d fn= kn*delta*n + 2.*gn*sqrt(0.5*me*kn)*vn;
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
		// double ftd = mu_d*fn.norm();
		// Tangential force
		gt = 0.2*2.*gn*sqrt(0.5*me*kn);
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
		pi->Fc += fn+ft;
		pj->Fc -= fn+ft;

		pi->Tc += (pi->R-0.5*delta)*arm;
		pj->Tc += (pj->R-0.5*delta)*arm;
	}
}

inline double DEM::EffectiveValue(double ai, double aj)
{
	double a = 0.;
	if (ai>1.0e-12 && aj>1.0e-12)	a = 2.*ai*aj/(ai+aj);
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
			if (Flag[ic][jc][kc]==-1)
			{
				Flag[ic][jc][kc] = p0->ID;
			}
			else
			{
				size_t q = Flag[ic][jc][kc];
				size_t p = p0->ID;
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
}

inline void DEM::FindContact()
{
	CMap.clear();
	for (size_t p=6; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];

		for (size_t d=0; d<3; ++d)
		{
			if (Periodic[d])
			{
				if (p0->X(d)+p0->R>DomSize[d])
				{
					p0->Xmir = p0->X;
					p0->Xmir(d) = p0->X(d)-DomSize[d]-1;
				}
				else if (p0->X(d)-p0->R<0.)			
				{
					p0->Xmir = p0->X;
					p0->Xmir(d) = p0->X(d)+DomSize[d]+1;
				}
			}
		}

		for (int i=p0->MinX; i<=p0->MaxX; ++i)
		for (int j=p0->MinY; j<=p0->MaxY; ++j)
		for (int k=p0->MinZ; k<=p0->MaxZ; ++k)
		{
			int ic = (i+Nx+1)%(Nx+1);
			int jc = (j+Ny+1)%(Ny+1);
			int kc = (k+Nz+1)%(Nz+1);
			
			int ind = Flag[ic][jc][kc];
			if (ind>5)
			{
				if (!Lp[ind]->fixed)	Flag[ic][jc][kc] = -1;
			}
		}
	}

	for (size_t p=6; p<Lp.size(); ++p)
	{
		DEM_PARTICLE* p0 = Lp[p];
		if (!p0->fixed)	UpdateFlag(p0);
	}


// used to double check contact
/*		for (size_t p=6; p<Lp.size(); ++p)
		for (size_t q=p+1; q<Lp.size(); ++q)
		{
			DEM_PARTICLE* pi = Lp[p];
			DEM_PARTICLE* pj = Lp[q];

			// Normal direction (pj pinnts to pi)
			Vector3d n = pi->X-pj->X;
			// Overlapping distance
			double delta = pi->R+pj->R-n.norm();

			Vector3d xi = pi->X;
			if (pi->X(0)+pi->R>Nx)	xi(0) -= Nx+1;
			if (pi->X(0)+pi->R<0)	xi(0) += Nx+1;

			Vector3d n1 = xi-pj->X;
			double delta1 = pi->R+pj->R-n1.norm();

			Vector3d xj = pj->X;
			if (pj->X(0)+pj->R>Nx)	xj(0) -= Nx+1;
			if (pj->X(0)+pj->R<0)	xj(0) += Nx+1;

			Vector3d n2 = pi->X-xj;
			double delta2 = pi->R+pj->R-n2.norm();

			if (delta>0. || delta1>0. || delta2>0.)
			{
				vector< size_t > par = {min(p,q), max(p,q)};
				if (find(Lc.begin(),Lc.end(),par)!=Lc.end())
				{

				}
				else
				{
					cout << "cannot find contact for par" << endl;
					abort();
				}
			}
		}

		for (size_t p=6; p<Lp.size(); ++p)
		{
			DEM_PARTICLE* p0 = Lp[p];
			if (p0->X(2)-p0->R<0.)
			{
				vector< size_t > par = {4,p};
				if (find(Lc.begin(),Lc.end(),par)!=Lc.end())
				{

				}
				else
				{
					cout << "cannot find contact for p and bot" << endl;
					abort();
				}
			}
			if (p0->X(1)-p0->R<0.)
			{
				vector< size_t > par = {2,p};
				if (find(Lc.begin(),Lc.end(),par)!=Lc.end())
				{

				}
				else
				{
					cout << "cannot find contact for p and ymin" << endl;
					abort();
				}
			}
			if (p0->X(1)+p0->R>Ny)
			{
				vector< size_t > par = {3,p};
				if (find(Lc.begin(),Lc.end(),par)!=Lc.end())
				{

				}
				else
				{
					cout << "cannot find contact for p and ymax" << endl;
					abort();
				}
			}
		}*/
}

inline void DEM::Contact()
{
    if (Lc.size()>0.)
    {
    	unordered_map<size_t, Vector3d> fmap;
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

inline void DEM::Solve(int tt, int ts, double dt)
{
	Dt = dt;
	for (int t=0; t<tt; ++t)
	{
		if (t%ts == 0)
		{
			cout << "Time Step = " << t << endl;
			WriteFileH5(t);
		}
		// clock_t t_start = std::clock();
		FindContact();
		// clock_t t_end = std::clock();
		// cout << "FindContact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		// t_start = std::clock();
		Contact();
		// t_end = std::clock();
		// cout << "Contact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (t==0)
		{
			for (size_t p=0; p<Lp.size(); ++p)
			{
				Lp[p]->Avb = (Lp[p]->Fh + Lp[p]->Fc + Lp[p]->Fex)/Lp[p]->M + Lp[p]->G;
				Lp[p]->Awb = Lp[p]->I.asDiagonal().inverse()*((Lp[p]->Th + Lp[p]->Tc + Lp[p]->Tex));
			}
		}
		// t_start = std::clock();
		Move();
		// t_end = std::clock();
		// cout << "Move time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		// t_start = std::clock();
		ZeroForceTorque();
		// t_end = std::clock();
		// cout << "ZeroForceTorque time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	}
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
	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));
    DataSet	*dataset_agv	= new DataSet(file.createDataSet("AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector));

    DataSet	*dataset_poi	= new DataSet(file.createDataSet("Points", PredType::NATIVE_DOUBLE, *space_points));
    DataSet	*dataset_fac	= new DataSet(file.createDataSet("Faces", PredType::NATIVE_INT, *space_faces));
    DataSet	*dataset_fv		= new DataSet(file.createDataSet("FaceVelocity", PredType::NATIVE_DOUBLE, *space_fscalar));

	dataset_r->write(r_h5, PredType::NATIVE_DOUBLE);
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
	delete dataset_tag;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_agv;

	delete dataset_poi;
	delete dataset_fac;
	delete dataset_fv;

	delete r_h5;
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
        ostringstream info;
        info << "particleInfo_" << p << ".res";
        ofstream log(info.str().c_str(), ios_base::out | ios_base::app);
        if (n==0)   log << "\"Time\"     \"X\"     \"Y\"     \"Z\"     \"U\"     \"V\"     \"W\"      \"Wx\"     \"Wy\"     \"Wz\"      \"Fhx\"     \"Fhy\"     \"Fhz\"      \"Thx\"     \"Thy\"     \"Thz\"     \"Fcx\"     \"Fcy\"     \"Fcz\"      \"Tcx\"     \"Tcy\"     \"Tcz\"\n";
        else;
        log << setprecision(9) << fixed << (double)n << "     " << Lp[p]->X(0) << "     " << Lp[p]->X(1) << "     " << Lp[p]->X(2) << "     " << Lp[p]->V(0) << "     " << Lp[p]->V(1) << "     " << Lp[p]->V(2) << "     " << Lp[p]->W(0) << "     " << Lp[p]->W(1) << "     " << Lp[p]->W(2) << "     " << Lp[p]->Fh(0) << "     " << Lp[p]->Fh(1) << "     " << Lp[p]->Fh(2)<< "     " << Lp[p]->Th(0) << "     " << Lp[p]->Th(1) << "     " << Lp[p]->Th(2) << "     " << Lp[p]->Fc(0) << "     " << Lp[p]->Fc(1) << "     " << Lp[p]->Fc(2)<< "     " << Lp[p]->Tc(0) << "     " << Lp[p]->Tc(1)<< "     " << Lp[p]->Tc(2) << "\n";
    }
}