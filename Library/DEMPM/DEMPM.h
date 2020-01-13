#include "../HEADER.h"
#include <DEM.h>
#include <MPM.h>

class DEMPM
{
public:
	DEMPM();
	~DEMPM();
	DEMPM(size_t nx, size_t ny, size_t nz, Vector3d dx, size_t ntype, string cmtype, string dmtype, double cr);
	void Init();
	void HertzContactPara(DEM_PARTICLE* pi, MPM_PARTICLE* pj, double delta, double& kn, double& gn);
	void Contact2P(DEM_PARTICLE* pi, MPM_PARTICLE* pj);
	// void ParticleToNodeWithDEM();
	void NodeToParticleWithDEM();
	void DEMtoNodeLocal(DEM_PARTICLE* p0);
	void DEMtoNode(bool show);
	void FindDEMContact();
	void ContactDEMP();
	void Solve(int tt, int ts);
	void WriteFileH5(int n);

	MPM*							DomMPM;												// Domain of MPM
	DEM*							DomDEM;												// Domain of DEM
	size_t 							Nnode;												// Total number of grids
    int 							Nx;													// Domain size
    int 							Ny;
    int 							Nz;	

	size_t 							Nproc;												// Number of processors which used
    int 							D;													// Demension
    vector < vector< size_t > >     Lc;                                                 // List of potential contacted paricles' ID
    vector <size_t>					LAn;												// List of actived nodes
};

inline DEMPM::DEMPM(size_t nx, size_t ny, size_t nz, Vector3d dx, size_t ntype, string cmtype, string dmtype, double cr)
{
	DomDEM = new DEM(nx, ny, nz, cmtype, dmtype, cr);
	DomMPM = new MPM(ntype, nx, ny, nz, dx);

	Nx = nx;
	Ny = ny;
	Nz = nz;
	Nnode = (Nx+1)*(Ny+1)*(Nz+1);
	D = DomMPM->D;
}

inline void DEMPM::Init()
{
	// DomDEM->Init();
	// DomMPM->Init();
	cout << "================ Start init.  ================" << endl;
	DomMPM->Lp.resize(0);
	DomMPM->Ln.resize(0);
	// Add basic nodes
	for (size_t n=0; n<Nnode; ++n)
	{
    	size_t i, j, k;
    	DomMPM->FindIndex(n, i, j, k);
		Vector3d x (i, j, k);
		DomMPM->Ln.push_back(new MPM_NODE(0,x));
    	DomMPM->Ln[DomMPM->Ln.size()-1]->ID = DomMPM->Ln.size()-1;
    	// set DPs_proc for coupling
    	DomMPM->Ln[DomMPM->Ln.size()-1]->DPs_proc.resize(DomDEM->Nproc);
	}

	DomDEM->Dt = 1.;
	Vector3d x0 (0., 0., 0.);
	DomDEM->Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
	DomDEM->Lp[0]->ID = 0;
	DomDEM->Lp[0]->Type = -1;
	DomDEM->Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
	DomDEM->Lp[1]->ID = 1;
	DomDEM->Lp[1]->Type = -1;
	DomDEM->Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
	DomDEM->Lp[2]->ID = 2;
	DomDEM->Lp[2]->Type = -1;
	DomDEM->Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
	DomDEM->Lp[3]->ID = 3;
	DomDEM->Lp[3]->Type = -1;
	DomDEM->Lp.push_back(new DEM_PARTICLE(-5, x0, 0.));
	DomDEM->Lp[4]->ID = 4;
	DomDEM->Lp[4]->Type = -1;
	DomDEM->Lp.push_back(new DEM_PARTICLE(-6, x0, 0.));
	DomDEM->Lp[5]->ID = 5;
	DomDEM->Lp[5]->Type = -1;
	cout << "=============== Finish init.  ================" << endl;
}

void DEMPM::HertzContactPara(DEM_PARTICLE* pi, MPM_PARTICLE* pj, double delta, double& kn, double& gn)
{
	double re 	= DomDEM->EffectiveValue(pi->R, pj->R);
	double me 	= DomDEM->EffectiveValue(pi->M, pj->M);
	double ee 	= 1./((1.-pi->Poisson*pi->Poisson)/pi->Young + (1.-pj->Poisson*pj->Poisson)/pj->Young);

	double sn 	= 2.*ee*sqrt(re*delta);
	kn 	= 1.3333333333333*ee*sqrt(re*delta);
	gn 	= -1.825741858351*DomDEM->Beta*sqrt(sn*me);
}

// Contact force model for DEM and MPM particles (only consider normal contact force for now)
inline void DEMPM::Contact2P(DEM_PARTICLE* pi, MPM_PARTICLE* pj)
{
	Vector3d Xi = pi->X;
	Vector3d Xj = pj->X;

	// Normal direction (pj pinnts to pi)
	Vector3d n = Xi-Xj;
	// Overlapping distance
	double delta = pi->R+pj->R-n.norm();

	if (delta>0.)
	{
		// if (delta>0.01*(pi->R+pj->R))
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
		n.normalize();
		double kn, gn;
		HertzContactPara(pi, pj, delta, kn, gn);
 		// Relative velocity in normal direction
		Vector3d vn = (pj->V-pi->V).dot(n)*n;
		// Normal contact force
		Vector3d fn= kn*delta*n + gn*vn;
		// cout << fn.transpose() << endl;
		// cout << n.transpose() << endl;
		// cout << "kn " << kn << endl;
		// cout << "delta " << delta << endl;
		// cout << "gn " << gn << endl;
		// cout << "vb " << vn.transpose() << endl;
		// abort();
        for (int d=0; d<D; ++d)
        {
            #pragma omp atomic
        	pi->Fc(d) += fn(d);
            #pragma omp atomic
            pj->Fc(d) -= fn(d);
        }		
	}
}

inline void DEMPM::DEMtoNodeLocal(DEM_PARTICLE* p0)
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
			int n = ic+jc*DomMPM->Ncy+kc*DomMPM->Ncz;
			size_t p = p0->ID;
			// if (DomMPM->Ln[n]->DPs.size()>0)
			// {
			// 	/*find DEM contact here*/
			// }
			DomMPM->Ln[n]->DPs.push_back(p);
			// cout << "n= " << n << endl;
			// cout << "p= " << p << endl;
			// cout << DomMPM->Ln[n]->DPs.size() << endl; 
			LAn.push_back(n);
		}
	}
}

void DEMPM::DEMtoNode(bool show)
{
	auto t_start = std::chrono::system_clock::now();
	// seperate Lc for openMP
	vector<vector <size_t>> lans;
	lans.resize(DomDEM->Nproc);
	#pragma omp parallel for schedule(static) num_threads(DomDEM->Nproc)
	for (size_t p=2*D; p<DomDEM->Lp.size(); ++p)
	{
		// cout << "p= " << p << endl;
		DEM_PARTICLE* p0 = DomDEM->Lp[p];
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
				// cout << "dis" << endl;
				int ic = (i+Nx+1)%(Nx+1);
				int jc = (j+Ny+1)%(Ny+1);
				int kc = (k+Nz+1)%(Nz+1);
				int n = ic+jc*DomMPM->Ncy+kc*DomMPM->Ncz;
				size_t p = p0->ID;

				auto proc_id = omp_get_thread_num();
				DomMPM->Ln[n]->DPs_proc[proc_id].push_back(p);
				lans[proc_id].push_back(n);
			}
		}
	}
	auto t_end = std::chrono::system_clock::now();
	if (show) 	cout << "000000000000000000= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
	// remove repeated elements for Lc
	#pragma omp parallel for schedule(static) num_threads(DomDEM->Nproc)
	for (size_t n=0; n<DomDEM->Nproc; ++n)
	{
		sort( lans[n].begin(), lans[n].end() );
		lans[n].erase(unique(lans[n].begin(), lans[n].end()), lans[n].end());
	}
	t_end = std::chrono::system_clock::now();
	if (show)	cout << "+++++++++++++++++++++= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	t_start = std::chrono::system_clock::now();
	for (size_t n=0; n<DomDEM->Nproc; ++n)
	{
		LAn.insert( LAn.end(), lans[n].begin(), lans[n].end() );
	}
	sort( LAn.begin(), LAn.end() );
	LAn.erase(unique(LAn.begin(), LAn.end()), LAn.end());
	t_end = std::chrono::system_clock::now();
	if (show)	cout << "//////////////////////= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	// cout << "2222222222" << endl;
	if (show)	cout << "LAn= " << LAn.size() << endl;
	// t_start = std::chrono::system_clock::now();
	// // combine DPs for collision
	// #pragma omp parallel for schedule(static) num_threads(DomDEM->Nproc+DomMPM->Nproc)
	// for (size_t n=0; n<LAn.size(); ++n)
	// {
	// 	size_t id = LAn[n];
	// 	DomMPM->Ln[id]->CombineDPs();
 //    }
 //    t_end = std::chrono::system_clock::now();
 //    if (show)	cout << "------------------= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
}

void DEMPM::FindDEMContact()
{
	// seperate Lc for openMP
	vector<vector<vector <size_t>>> lcs;
	lcs.resize(DomDEM->Nproc+DomMPM->Nproc);
	// combine DPs
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<LAn.size(); ++n)
	{
		size_t id = LAn[n];
		DomMPM->Ln[id]->CombineDPs();
		if (DomMPM->Ln[id]->DPs.size()>1)
		{
			for (size_t i=0; i<DomMPM->Ln[id]->DPs.size(); ++i)
			for (size_t j=i+1; j<DomMPM->Ln[id]->DPs.size(); ++j)
			{
				size_t p = DomMPM->Ln[id]->DPs[i];
				size_t q = DomMPM->Ln[id]->DPs[j];
				auto proc_id = omp_get_thread_num();
				lcs[proc_id].push_back({min(p,q), max(p,q)});
				if (p==q)
				{
					#pragma omp critical
					{
						cout << "p=q" << endl;
						for (size_t k=0; k<DomMPM->Ln[id]->DPs.size(); ++k)
						{
							cout << "dps: " << DomMPM->Ln[id]->DPs[k] << endl;
						}
						// cout << DomMPM->Ln[id]->DPs.size() << endl;
						// cout << p << endl;
						abort();
					}
				}
				// #pragma omp critical
				// {
				// 	cout << p << endl;
				// 	cout << q << endl;
				// 	// abort();
				// }
			}
		}
    }
    // cout << "sssssssssssssssssss" << endl;
    // abort();
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t n=0; n<Nproc; ++n)
	{
		sort( lcs[n].begin(), lcs[n].end() );
		lcs[n].erase(unique(lcs[n].begin(), lcs[n].end()), lcs[n].end());
	}
	for (size_t n=0; n<Nproc; ++n)
	{
		DomDEM->Lc.insert( DomDEM->Lc.end(), lcs[n].begin(), lcs[n].end() );
	}
	sort( DomDEM->Lc.begin(), DomDEM->Lc.end() );
	DomDEM->Lc.erase(unique(DomDEM->Lc.begin(), DomDEM->Lc.end()), DomDEM->Lc.end());
	// cout << "bbbbbbbbbbbbbbbbbbbb" << endl;
}

// void DEMPM::DEMtoNode()
// {
// 	for (size_t p=6; p<DomDEM->Lp.size(); ++p)
// 	{
// 		DEMtoNodeLocal(DomDEM->Lp[p]);
// 	}
// }

void DEMPM::NodeToParticleWithDEM()
{
	auto t_start = std::chrono::system_clock::now();
	// seperate Lc for openMP
	vector<vector<vector <size_t>>> lcs;
	lcs.resize(DomMPM->Nproc);
	#pragma omp parallel for schedule(static) num_threads(DomMPM->Nproc)
	for (size_t p=0; p<DomMPM->Lp.size(); ++p)
	{
		// Reset position increasement of this particle
		DomMPM->Lp[p]->DeltaX 	= Vector3d::Zero();
		DomMPM->Lp[p]->StressSmooth.setZero();
		if (!DomMPM->Lp[p]->FixV)
		{
			for (size_t l=0; l<DomMPM->Lp[p]->Lni.size(); ++l)
			{
				size_t 	id = DomMPM->Lp[p]->Lni[l];
				double 	n  = DomMPM->Lp[p]->LnN[l];
				// Update velocity of this particle
				Vector3d an = DomMPM->Ln[id]->F/DomMPM->Ln[id]->M;
				DomMPM->Lp[p]->V += n*an*DomMPM->Dt;
				DomMPM->Lp[p]->X += n*DomMPM->Ln[id]->V*DomMPM->Dt;
				DomMPM->Lp[p]->StressSmooth += n*DomMPM->Ln[id]->Stress;
			}
		}
		else
		{
			DomMPM->Lp[p]->V = DomMPM->Lp[p]->Vf;
			DomMPM->Lp[p]->X += DomMPM->Lp[p]->V*DomMPM->Dt;
		}
		// Velocity gradient tensor
		DomMPM->CalVGradLocal(p);
		// Update deformation tensor
		DomMPM->Lp[p]->F = (Matrix3d::Identity() + DomMPM->Lp[p]->L*DomMPM->Dt)*DomMPM->Lp[p]->F;
		// Update particle length
		if (DomMPM->Lp[p]->Type==0)	DomMPM->CalPSizeR(p);
		// CalPSizeCP(p);
		// Update volume of particles
		DomMPM->Lp[p]->Vol 	= DomMPM->Lp[p]->F.determinant()*DomMPM->Lp[p]->Vol0;
		// Update strain
		Matrix3d de = 0.5*DomMPM->Dt*(DomMPM->Lp[p]->L + DomMPM->Lp[p]->L.transpose());
		// Update stress
		Matrix3d w = 0.5*DomMPM->Dt*((DomMPM->Lp[p]->L - DomMPM->Lp[p]->L.transpose()));
		DomMPM->Lp[p]->Stress += w*DomMPM->Lp[p]->Stress - DomMPM->Lp[p]->Stress*w.transpose();
		if (DomMPM->Lp[p]->Type==0)			DomMPM->Lp[p]->Elastic(de);
		else if (DomMPM->Lp[p]->Type==1)
		{
			DomMPM->Lp[p]->EOSMonaghan(DomMPM->Cs);
			DomMPM->Lp[p]->Newtonian(de);
		}
		else if (DomMPM->Lp[p]->Type==2)	DomMPM->Lp[p]->MohrCoulomb(de);
		else if (DomMPM->Lp[p]->Type==3)	DomMPM->Lp[p]->DruckerPrager(de);
		// Reset hydro force and contact force
		DomMPM->Lp[p]->Fh.setZero();
		DomMPM->Lp[p]->Fc.setZero();
		// cout << "start check mpm particle range" << endl;
		// find contact between DEM and MPM
		int maxx = (int) (DomMPM->Lp[p]->X(0)+DomMPM->Lp[p]->R+1.);
		int maxy = (int) (DomMPM->Lp[p]->X(1)+DomMPM->Lp[p]->R+1.);
		int maxz = (int) (DomMPM->Lp[p]->X(2)+DomMPM->Lp[p]->R+1.);

		int minx = (int) (DomMPM->Lp[p]->X(0)-DomMPM->Lp[p]->R-1.);
		int miny = (int) (DomMPM->Lp[p]->X(1)-DomMPM->Lp[p]->R-1.);
		int minz = (int) (DomMPM->Lp[p]->X(2)-DomMPM->Lp[p]->R-1.);
		// for surounding box
		for (int i=minx; i<=maxx; ++i)
		for (int j=miny; j<=maxy; ++j)
		for (int k=minz; k<=maxz; ++k)
		{
			Vector3d ind (i,j,k);

			int n = i+j*DomMPM->Ncy+k*DomMPM->Ncz;
			if (DomMPM->Ln[n]->DPs.size()>0)
			{
				// distance between node to MPM particles
				double dis = (ind-DomMPM->Lp[p]->X).norm()-DomMPM->Lp[p]->R-0.87;
				if (dis<0.)
				{
					auto proc_id = omp_get_thread_num();
					// potential contact between MPM and DEM
					for (size_t i=0; i<DomMPM->Ln[n]->DPs.size(); ++i)
					{
						size_t q = DomMPM->Ln[n]->DPs[i];
						lcs[proc_id].push_back({q,p});
					}
					// #pragma omp critical
					// {
					// 	// potential contact between MPM and DEM
					// 	for (size_t i=0; i<DomMPM->Ln[n]->DPs.size(); ++i)
					// 	{
					// 		size_t q = DomMPM->Ln[n]->DPs[i];
					// 		Lc.push_back({q,p});
					// 	}
					// }
				}
			}
		}
	}
	auto t_end = std::chrono::system_clock::now();
	// cout << "111111111111= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
	// remove repeated elements for Lc
	#pragma omp parallel for schedule(static) num_threads(DomMPM->Nproc)
	for (size_t n=0; n<DomMPM->Nproc; ++n)
	{
		sort( lcs[n].begin(), lcs[n].end() );
		lcs[n].erase(unique(lcs[n].begin(), lcs[n].end()), lcs[n].end());
	}
	t_end = std::chrono::system_clock::now();
	// cout << "2222222222= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
	t_start = std::chrono::system_clock::now();
	for (size_t n=0; n<DomMPM->Nproc; ++n)
	{
		Lc.insert( Lc.end(), lcs[n].begin(), lcs[n].end() );
	}
	sort( Lc.begin(), Lc.end() );
	Lc.erase(unique(Lc.begin(), Lc.end()), Lc.end());	
	t_end = std::chrono::system_clock::now();	
	// cout << "3333333333= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
}

void DEMPM::ContactDEMP()
{
	// make Lc unique
	// sort( LAn.begin(), LAn.end() );
	// LAn.erase(unique(LAn.begin(), LAn.end()), LAn.end());
	sort( Lc.begin(), Lc.end() );
	Lc.erase(unique(Lc.begin(), Lc.end()), Lc.end());
    if (Lc.size()>0.)
    {
    	// #pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t l=0; l<Lc.size(); ++l)
		{
			int i = Lc[l][0];
			int j = Lc[l][1];
			Contact2P(DomDEM->Lp[i], DomMPM->Lp[j]);
		}
    }
}

void DEMPM::Solve(int tt, int ts)
{
	for (int t=0; t<tt; ++t)
	{
		bool show = false;
		if (t%ts==0)	show = true;
		if (show) 	cout << "Time Step = " << t << endl;
		if (t%ts == 0)
		// if (t> 2000)
		{
			WriteFileH5(t);
		}
		if (show)	cout << "start contact" << endl;
		ContactDEMP();
		if (show)	cout << "finish contact" << endl;
		Lc.clear();
		if (show)	cout << "start move" << endl;
		DomDEM->Move();
		if (show)	cout << "finish move" << endl;
		DomDEM->ZeroForceTorque(true, true);
		for (size_t n=0; n<LAn.size(); ++n)
		{
			size_t id = LAn[n];
			DomMPM->Ln[id]->ResetwithDEM(DomDEM->Nproc);
	    }
	    LAn.resize(0);
	    if (show)	cout << "finish ZeroForceTorque" << endl;
		// cout << LAn.size() << endl;
		// cout << "11111111111111111" << endl;
		auto t_start = std::chrono::system_clock::now();
		DomMPM->ParticleToNode();
		auto t_end = std::chrono::system_clock::now();
		if (show)	cout << "ParticleToNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		DomMPM->CalVOnNode();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "CalVOnNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		DEMtoNode(show);
		t_end = std::chrono::system_clock::now();
		// if (show)	cout << "finish DEMtoNode" << endl;
		if (show)	cout << "DEMtoNode= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		FindDEMContact();
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "FindDEMContact= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		DomDEM->Contact(false, 0);
		int nc = DomDEM->Lc.size();
		DomDEM->Lc.clear();
		t_end = std::chrono::system_clock::now();
		if (show)   cout << "nc= " << nc << endl;
		if (show)	cout << "DEM Contact time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		t_start = std::chrono::system_clock::now();
		NodeToParticleWithDEM();
		// if (Lc.size()>0) 	cout << "Lc.size()= " << Lc.size() << endl;
		t_end = std::chrono::system_clock::now();
		if (show)	cout << "NodeToParticle= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
		if (show) 	cout << "===========================" << endl;
	}
}

inline void DEMPM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = "DEMPM_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar_mpm[1] = {DomMPM->Lp.size()};			//create data space.
	hsize_t	dims_vector_mpm[1] = {3*DomMPM->Lp.size()};			//create data space.
	hsize_t	dims_tensor_mpm[1] = {6*DomMPM->Lp.size()};

	int rank_scalar_mpm = sizeof(dims_scalar_mpm) / sizeof(hsize_t);
	int rank_vector_mpm = sizeof(dims_vector_mpm) / sizeof(hsize_t);
	int rank_tensor_mpm = sizeof(dims_tensor_mpm) / sizeof(hsize_t);

	DataSpace	*space_scalar_mpm = new DataSpace(rank_scalar_mpm, dims_scalar_mpm);
	DataSpace	*space_vector_mpm = new DataSpace(rank_vector_mpm, dims_vector_mpm);
	DataSpace	*space_tensor_mpm = new DataSpace(rank_tensor_mpm, dims_tensor_mpm);

	double* tag_h5_mpm 	= new double[  DomMPM->Lp.size()];
	double* m_h5_mpm 	= new double[  DomMPM->Lp.size()];
	double* you_h5_mpm 	= new double[  DomMPM->Lp.size()];
	double* poi_h5_mpm 	= new double[  DomMPM->Lp.size()];
	double* pos_h5_mpm 	= new double[3*DomMPM->Lp.size()];
	double* vel_h5_mpm 	= new double[3*DomMPM->Lp.size()];
	double* s_h5_mpm	= new double[6*DomMPM->Lp.size()];

	for (size_t i=0; i<DomMPM->Lp.size(); ++i)
	{
        tag_h5_mpm[  i  ] 	= DomMPM->Lp[i]->Tag;
        m_h5_mpm  [  i  ] 	= DomMPM->Lp[i]->M;
        you_h5_mpm[  i  ] 	= DomMPM->Lp[i]->Young;
        poi_h5_mpm[  i  ] 	= DomMPM->Lp[i]->Poisson;
		pos_h5_mpm[3*i  ] 	= DomMPM->Lp[i]->X(0);
		pos_h5_mpm[3*i+1] 	= DomMPM->Lp[i]->X(1);
		pos_h5_mpm[3*i+2] 	= DomMPM->Lp[i]->X(2);
		vel_h5_mpm[3*i  ] 	= DomMPM->Lp[i]->V(0);
		vel_h5_mpm[3*i+1] 	= DomMPM->Lp[i]->V(1);
		vel_h5_mpm[3*i+2] 	= DomMPM->Lp[i]->V(2);
		s_h5_mpm  [6*i  ] 	= DomMPM->Lp[i]->StressSmooth(0,0);
		s_h5_mpm  [6*i+1] 	= DomMPM->Lp[i]->StressSmooth(0,1);
		s_h5_mpm  [6*i+2] 	= DomMPM->Lp[i]->StressSmooth(0,2);
		s_h5_mpm  [6*i+3] 	= DomMPM->Lp[i]->StressSmooth(1,1);
		s_h5_mpm  [6*i+4] 	= DomMPM->Lp[i]->StressSmooth(1,2);
		s_h5_mpm  [6*i+5] 	= DomMPM->Lp[i]->StressSmooth(2,2);
	}

	DataSet	*dataset_tag_mpm	= new DataSet(file.createDataSet("MPM_Tag", PredType::NATIVE_DOUBLE, *space_scalar_mpm));
	DataSet	*dataset_m_mpm		= new DataSet(file.createDataSet("MPM_Mass", PredType::NATIVE_DOUBLE, *space_scalar_mpm));
	DataSet	*dataset_you_mpm	= new DataSet(file.createDataSet("MPM_Young", PredType::NATIVE_DOUBLE, *space_scalar_mpm));
	DataSet	*dataset_poi_mpm	= new DataSet(file.createDataSet("MPM_Poisson", PredType::NATIVE_DOUBLE, *space_scalar_mpm));
    DataSet	*dataset_pos_mpm	= new DataSet(file.createDataSet("MPM_Position", PredType::NATIVE_DOUBLE, *space_vector_mpm));
    DataSet	*dataset_vel_mpm	= new DataSet(file.createDataSet("MPM_Velocity", PredType::NATIVE_DOUBLE, *space_vector_mpm));
    DataSet	*dataset_s_mpm		= new DataSet(file.createDataSet("MPM_Stress", PredType::NATIVE_DOUBLE, *space_tensor_mpm));

	dataset_tag_mpm->write(tag_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_m_mpm->write(m_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_you_mpm->write(you_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_poi_mpm->write(poi_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_pos_mpm->write(pos_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_vel_mpm->write(vel_h5_mpm, PredType::NATIVE_DOUBLE);
	dataset_s_mpm->write(s_h5_mpm, PredType::NATIVE_DOUBLE);

	delete space_scalar_mpm;
	delete space_vector_mpm;
	delete space_tensor_mpm;
	delete dataset_tag_mpm;
	delete dataset_m_mpm;
	delete dataset_you_mpm;
	delete dataset_poi_mpm;
	delete dataset_pos_mpm;
	delete dataset_vel_mpm;
	delete dataset_s_mpm;

	delete tag_h5_mpm;
	delete m_h5_mpm;
	delete you_h5_mpm;
	delete poi_h5_mpm;
	delete pos_h5_mpm;
	delete vel_h5_mpm;
	delete s_h5_mpm;

	// DEM datas 
	size_t np_dem = DomDEM->Lp.size()-2*D;
	hsize_t	dims_scalar_dem[1] = {np_dem};			//create data space.
	hsize_t	dims_vector_dem[1] = {3*np_dem};			//create data space.

	int rank_scalar_dem = sizeof(dims_scalar_dem) / sizeof(hsize_t);
	int rank_vector_dem = sizeof(dims_vector_dem) / sizeof(hsize_t);

	DataSpace	*space_scalar_dem = new DataSpace(rank_scalar_dem, dims_scalar_dem);
	DataSpace	*space_vector_dem = new DataSpace(rank_vector_dem, dims_vector_dem);

	double* r_h5_dem 	= new double[  np_dem];
	double* rho_h5_dem 	= new double[  np_dem];
	double* tag_h5_dem 	= new double[  np_dem];
	double* pos_h5_dem 	= new double[3*np_dem];
	double* vel_h5_dem 	= new double[3*np_dem];
	double* agv_h5_dem	= new double[3*np_dem];
	double* fh_h5_dem 	= new double[3*np_dem];

	for (size_t i=0; i<np_dem; ++i)
	{
		size_t ip = i+2*D;
		Vector3d agv = DomDEM->Lp[ip]->Q._transformVector(DomDEM->Lp[ip]->W);
        r_h5_dem  [  i  ] 	= DomDEM->Lp[ip]->R;
        rho_h5_dem[  i  ] 	= DomDEM->Lp[ip]->Rho;
        tag_h5_dem[  i  ] 	= DomDEM->Lp[ip]->Tag;
		pos_h5_dem[3*i  ] 	= DomDEM->Lp[ip]->X(0);
		pos_h5_dem[3*i+1] 	= DomDEM->Lp[ip]->X(1);
		pos_h5_dem[3*i+2] 	= DomDEM->Lp[ip]->X(2);
		vel_h5_dem[3*i  ] 	= DomDEM->Lp[ip]->V(0);
		vel_h5_dem[3*i+1] 	= DomDEM->Lp[ip]->V(1);
		vel_h5_dem[3*i+2] 	= DomDEM->Lp[ip]->V(2);
		fh_h5_dem[3*i  ] 	= DomDEM->Lp[ip]->Fh(0);
		fh_h5_dem[3*i+1] 	= DomDEM->Lp[ip]->Fh(1);
		fh_h5_dem[3*i+2] 	= DomDEM->Lp[ip]->Fh(2);
		agv_h5_dem[3*i  ] 	= agv(0);
		agv_h5_dem[3*i+1] 	= agv(1);
		agv_h5_dem[3*i+2] 	= agv(2);
	}

	DataSet	*dataset_r_dem 		= new DataSet(file.createDataSet("DEM_Radius", PredType::NATIVE_DOUBLE, *space_scalar_dem));
	DataSet	*dataset_rho_dem	= new DataSet(file.createDataSet("DEM_Rho", PredType::NATIVE_DOUBLE, *space_scalar_dem));
	DataSet	*dataset_tag_dem	= new DataSet(file.createDataSet("DEM_Tag", PredType::NATIVE_DOUBLE, *space_scalar_dem));
    DataSet	*dataset_pos_dem	= new DataSet(file.createDataSet("DEM_Position", PredType::NATIVE_DOUBLE, *space_vector_dem));
    DataSet	*dataset_vel_dem	= new DataSet(file.createDataSet("DEM_Velocity", PredType::NATIVE_DOUBLE, *space_vector_dem));
    DataSet	*dataset_agv_dem	= new DataSet(file.createDataSet("DEM_AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector_dem));
    DataSet	*dataset_fh_dem		= new DataSet(file.createDataSet("DEM_HydroForce", PredType::NATIVE_DOUBLE, *space_vector_dem));

	dataset_r_dem->write(r_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_rho_dem->write(rho_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_tag_dem->write(tag_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_pos_dem->write(pos_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_vel_dem->write(vel_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_agv_dem->write(agv_h5_dem, PredType::NATIVE_DOUBLE);
	dataset_fh_dem->write(fh_h5_dem, PredType::NATIVE_DOUBLE);

	delete space_scalar_dem;
	delete space_vector_dem;

	delete dataset_r_dem;
	delete dataset_rho_dem;
	delete dataset_tag_dem;
	delete dataset_pos_dem;
	delete dataset_vel_dem;
	delete dataset_agv_dem;
	delete dataset_fh_dem;

	delete r_h5_dem;
	delete rho_h5_dem;
	delete tag_h5_dem;
	delete pos_h5_dem;
	delete vel_h5_dem;
	delete agv_h5_dem;
	delete fh_h5_dem;

	file.close();

	string file_name_xmf = "DEMPM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << DomMPM->Lp.size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << DomMPM->Lp.size() << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/MPM_Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DomMPM->Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/MPM_Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DomMPM->Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/MPM_Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DomMPM->Lp.size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/MPM_Stress\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << np_dem << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << np_dem << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/DEM_Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_Rho \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_AngularVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"HydroForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << np_dem << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/DEM_HydroForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}