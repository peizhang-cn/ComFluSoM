#include "../HEADER.h"
#include <DEM.h>
#include <LBM.h>
#include <RWM.h>

class DELBM
{
public:
	DELBM();
	~DELBM();
	DELBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, string cmtype, string dmtype, double cr);
	void Init(double rho0, Vector3d initV);
    void SetRW(double dc, double dt);
    void InitBoundaryG();                                                               // Add boundaries as DEM particles and set G for them
	void InitG(DEM_PARTICLE* p0);
    void ResetG();
	void UpdateG();
    void UpdateGForIBB();
    void UpdateGForVAM();
	void UpdateLn(DEM_PARTICLE* p0);
    void Refill();
    void EqRefill();
    void NormRefill();
    void ApplyIBB(int method);                                                          // "0" for VIBB, "1" for Yu's double IBB
    void ApplyNEBB();
    void ApplyLNEBB();
    void ApplyVAM();
    void MoveRW();
    void MoveVAM();
    void UpdateBoxGlobal();
    void UpdateXbrForRWM();
    void CrossPBC(DEM_PARTICLE* p0);
    void ApplyPSM();
    double CalCrossPointQ(Vector3d& a, Vector3d& v, Vector3d& c, double r);
    double CalGdistance(int i, int j, int k, DEM_PARTICLE* p0);
    double CalGvolumeSphere(int i, int j, int k, DEM_PARTICLE* p0, int n);
	void SolveOneStepIBB(int method, int demNt);
    void SolveOneStepNEBB(int demNt);
    void SolveOneStepVAM(int demNt, bool save, int ct);
    void SolveOneStepPSM(int demNt);
    void SolveVAM(int tt, int savet, int demNt);
    void CalRhoVDELBM();
    void WriteFileH5(int n, int scale);


	LBM* 							DomLBM;												// Domain of LBM
	DEM*							DomDEM;												// Domain of DEM
    RWM*                            DomRWM;                                             // Domain of RWM

    bool                            Periodic[3];
    bool                            UseRW;
    int 							Nx;													// Domain size
    int 							Ny;
    int 							Nz;
    int                             DomSize[3];
    // int                             PeriodicIndex[3];

	int 							Nproc;												// Number of processors which used
    int 							D;

    double                          Hi;                                                 // Interpolation distance
    int                             BCmethod;                                           // Boundary condition methods, 0 for NEBB and 1 for IMB
    int                             Gmethod;                                            // Methods for calculate Gamma for IMB, 0 for distance based method and 1 for subgrid method

    vector< vector<int> >			Lr;                                                 // List of refilling nodes
};

inline DELBM::DELBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu, string cmtype, string dmtype, double cr)
{
	DomLBM = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, nu);
	DomDEM = new DEM(nx, ny, nz, cmtype, dmtype, cr);
    DomRWM = new RWM(nx, ny, nz, 0., 1.);
	Nx = nx;
	Ny = ny;
	Nz = nz;

    DomSize[0] = Nx;
    DomSize[1] = Ny;
    DomSize[2] = Nz;

    Periodic[0] = true;
    Periodic[1] = true;
    Periodic[2] = true;
    UseRW = false;
    D = DomLBM->D;

    Hi = 0.25;

    Nproc = 12;

    Vector3d x0 (0., 0., 0.);
    DomDEM->Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
    DomDEM->Lp[0]->ID = 0;
    DomDEM->Lp[0]->Type = -1;
    DomDEM->Lp[0]->X << 0., 0.5*Ny, 0.5*Nz;
    DomDEM->Lp[0]->Normal << 1., 0., 0.;
    DomDEM->Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
    DomDEM->Lp[1]->ID = 1;
    DomDEM->Lp[1]->Type = -1;
    DomDEM->Lp[1]->X << Nx, 0.5*Ny, 0.5*Nz;
    DomDEM->Lp[1]->Normal << -1., 0., 0.;

    DomDEM->Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
    DomDEM->Lp[2]->ID = 2;
    DomDEM->Lp[2]->Type = -1;
    DomDEM->Lp[2]->X << 0.5*Nx, 0., 0.5*Nz;
    DomDEM->Lp[2]->Normal << 0., 1., 0.;
    DomDEM->Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
    DomDEM->Lp[3]->ID = 3;
    DomDEM->Lp[3]->Type = -1;
    DomDEM->Lp[3]->X << 0.5*Nx, Ny, 0.5*Nz;
    DomDEM->Lp[3]->Normal << 0., -1., 0.;

    DomDEM->Lp.push_back(new DEM_PARTICLE(-5, x0, 0.));
    DomDEM->Lp[4]->ID = 4;
    DomDEM->Lp[4]->Type = -1;
    DomDEM->Lp[4]->X << 0.5*Nx, 0.5*Ny, 0.;
    DomDEM->Lp[4]->Normal << 0., 0., 1.;
    DomDEM->Lp.push_back(new DEM_PARTICLE(-6, x0, 0.));
    DomDEM->Lp[5]->ID = 5;
    DomDEM->Lp[5]->Type = -1;
    DomDEM->Lp[5]->X << 0.5*Nx, 0.5*Ny, Nz;
    DomDEM->Lp[5]->Normal << 0., 0., -1.;
}
// this function must be called before Init
inline void DELBM::SetRW(double dc, double dt)
{
    DomRWM->Dc = dc;
    DomRWM->Dt = dt;
    UseRW = true;
    // DomRWM->V_ptr = &DomLBM->V;
}

inline void DELBM::Init(double rho0, Vector3d initV)
{
	DomLBM->Init(rho0, initV);

    // Mark flag and add particles for the surounding box to handle walls and periodic BC
    // For x axis
    for (int j=0; j<=Ny; ++j)
    for (int k=0; k<=Nz; ++k)
    {
        DomLBM->G[0 ][j][k][0] = 0.;
        DomLBM->G[Nx][j][k][0] = 1.;
        DomLBM->Flag[0 ][j][k].push_back(0);
        DomLBM->Flag[Nx][j][k].push_back(1);
    }
    // For y axis
    for (int i=0; i<=Nx; ++i)
    for (int k=0; k<=Nz; ++k)
    {
        DomLBM->G[i][0 ][k][0] = 2.;
        DomLBM->G[i][Ny][k][0] = 3.;
        DomLBM->Flag[i][0 ][k].push_back(2);
        DomLBM->Flag[i][Ny][k].push_back(3);
    }       
    if (D>2)
    {
        // For z axis
        for (int i=0; i<=Nx; ++i)
        for (int j=0; j<=Ny; ++j)
        {
            DomLBM->G[i][j][0 ][0] = 4.;
            DomLBM->G[i][j][Nz][0] = 5.;
            DomLBM->Flag[i][j][0 ].push_back(4);
            DomLBM->Flag[i][j][Nz].push_back(5);
        }
    }

	for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        cout << "p= " << p << endl;
        InitG(DomDEM->Lp[p]);
    }

    DomLBM->Nproc = Nproc;
    DomRWM->Nproc = Nproc;
    DomRWM->Nproc = 1;

    Lr.resize(0);
    DomDEM->Lc.resize(0);

    if (UseRW)  
    {
        DomRWM->Init();
        DomRWM->V_ptr = DomLBM->V;
    }
    cout << "init finish" << endl;
}

// inline void DELBM::InitBoundaryG()
// {
//     // Mark flag and add particles for the surounding box to handle walls and periodic BC
//     // For x axis
//     for (int j=0; j<=Ny; ++j)
//     for (int k=0; k<=Nz; ++k)
//     {
//         DomLBM->G[0 ][j][k][0] = 0.;
//         DomLBM->G[Nx][j][k][0] = 1.;
//     }
//     Vector3d x0 (0., 0., 0.);
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-1, x0, 0.));
//     DomDEM->Lp[0]->ID = 0;
//     DomDEM->Lp[0]->Type = -1;
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-2, x0, 0.));
//     DomDEM->Lp[1]->ID = 1;
//     DomDEM->Lp[1]->Type = -1;
//     // For y axis
//     for (int i=0; i<=Nx; ++i)
//     for (int k=0; k<=Nz; ++k)
//     {
//         DomLBM->G[i][0 ][k][0] = 2.;
//         DomLBM->G[i][Ny][k][0] = 3.;
//     }       
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-3, x0, 0.));
//     DomDEM->Lp[2]->ID = 2;
//     DomDEM->Lp[2]->Type = -1;
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-4, x0, 0.));
//     DomDEM->Lp[3]->ID = 3;
//     DomDEM->Lp[3]->Type = -1;
//     // For z axis
//     for (int i=0; i<=Nx; ++i)
//     for (int j=0; j<=Ny; ++j)
//     {
//         DomLBM->G[i][j][0 ][0] = 4.;
//         DomLBM->G[i][j][Nz][0] = 5.;
//     }
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-5, x0, 0.));
//     DomDEM->Lp[4]->ID = 4;
//     DomDEM->Lp[4]->Type = -1;
//     DomDEM->Lp.push_back(new DEM_PARTICLE(-6, x0, 0.));
//     DomDEM->Lp[5]->ID = 5;
//     DomDEM->Lp[5]->Type = -1;
// }

inline void DELBM::InitG(DEM_PARTICLE* p0)
{
    int xMin, xMax, yMin, yMax, zMin, zMax;
    // Only init Gamma inside of the surrounding box for particles.
    xMin = p0->Min(0);
    xMax = p0->Max(0);

    yMin = p0->Min(1);
    yMax = p0->Max(1);

    zMin = p0->Min(2);
    zMax = p0->Max(2);
    // To aviod nagative index, need to rewrite for MPI.
    xMin = (xMin<0 ) ? 0 :xMin;
    xMax = (xMax>Nx) ? Nx:xMax;

    yMin = (yMin<0 ) ? 0 :yMin;
    yMax = (yMax>Ny) ? Ny:yMax;

    zMin = (zMin<0 ) ? 0 :zMin;
    zMax = (zMax>Nz) ? Nz:zMax;
    // For NEBB
    if (BCmethod==0)
    {
        for (int i=xMin; i<=xMax; ++i)
        for (int j=yMin; j<=yMax; ++j)
        for (int k=zMin; k<=zMax; ++k)
        {
            bool mark = false;
            Vector3d ind (i,j,k);
            double dis = (ind-p0->X).norm()-p0->R;
            double id = (double) p0->ID;
            double g = DomLBM->G[i][j][k][0];

            if (dis<0.)     mark = true;
            // First check for possible crossing, sqrt(3) < 1.75 (for 3d)
            else if (dis<1.75)
            {
                for (int q=0; q<DomLBM->Q; ++q)
                {
                    // Distance along normal direction
                    double disn = (ind+DomLBM->E[q]-p0->X).norm()-p0->R;
                    if (disn<0.)
                    {
                        mark = true;
                        break;
                    }
                }
                if (mark)
                {
                    // Reference point position
                    Vector3d xi = Hi*(ind-p0->X).normalized()+ind;
                    p0->Lb.push_back({i,j,k});
                    p0->Ld.push_back(dis);
                    p0->Li.push_back(xi);
                    // Mark boundary node as ID+0.5
                    id += 0.5;
                }
            }

            if (mark)
            {
                // If this node was fluid
                if (g==-1.)     DomLBM->G[i][j][k][0] = id;
                else
                {
                    // If this node was a boundary node and it will change to a solid node
                    if ((g-((int) g)) == 0.5 && (id-((int) id))==0)
                    {
                        DomLBM->G[i][j][k][0] = id;
                    }
                    size_t p = p0->ID;
                    size_t q = (int) DomLBM->G[i][j][k][0];
                    DomDEM->Lc.push_back({min(p,q), max(p,q)});
                }
            }
        }
    }
    // For Volume Average Methods
    else if (BCmethod==1)
    {
        cout << "Volume Average Methods" << endl;
        for (int i=xMin; i<=xMax; ++i)
        for (int j=yMin; j<=yMax; ++j)
        for (int k=zMin; k<=zMax; ++k)
        {
            // cout << i << " " << j << " " << k << endl;
            double tempG = 0;
            if (Gmethod==0)            tempG = CalGdistance(i,j,k,p0);
            else if (Gmethod==1)       tempG = CalGvolumeSphere(i,j,k,p0, 20);
            else
            {
                cout << "wrong Gmethod!" << endl;
                abort();
            }
            if (tempG>0)
            {
                // cout << "tempg" << endl;
                if (tempG<1.)
                {
                    p0->Lb.push_back({i,j,k});
                }
                if (DomLBM->G[i][j][k][0]==-1.)
                {
                    DomLBM->G[i][j][k][0] = p0->ID;
                    DomLBM->G[i][j][k][1] = tempG;
                }
                else
                {
                    DomLBM->G[i][j][k][2] = p0->ID;
                    DomLBM->G[i][j][k][3] = tempG;
                    size_t p0 = (int) DomLBM->G[i][j][k][0];
                    size_t p1 = DomDEM->Lp[p0]->ID;
                    DomDEM->Lc.push_back({min(p0,p1), max(p0,p1)});
                }
            }
        }
        for (int i=xMin; i<=xMax; ++i)
        for (int j=yMin; j<=yMax; ++j)
        for (int k=zMin; k<=zMax; ++k)
        {
            p0->Ln.push_back({i,j,k});
        }
    }
    else if (BCmethod==2)
    {
        for (int i=xMin; i<=xMax; ++i)
        for (int j=yMin; j<=yMax; ++j)
        for (int k=zMin; k<=zMax; ++k)
        {
            Vector3d ind (i,j,k);
            double dis = (ind-p0->X).norm()-p0->R;
            bool mark = false;
            double id = (double) p0->ID;
            if (dis<0.)     mark = true;
            // First check for possible crossing, sqrt(3) < 1.75 (for 3d)
            else if (dis<1.75)  
            {
                VectorXd deltaQ(DomLBM->Q);
                for (int q=0; q<DomLBM->Q; ++q)
                {
                    deltaQ(q) = CalCrossPointQ(ind, DomLBM->E[DomLBM->Op[q]], p0->X, p0->R);
                    if (deltaQ(q)>-1.)  mark = true;
                }

                if (mark)
                {
                    p0->Lb.push_back({i,j,k});
                    p0->Lq.push_back(deltaQ);
                    id = p0->ID+0.5;
                }
            } 
            if (mark)
            {
                double g = DomLBM->G[i][j][k][0];
                if (g==-1.)     DomLBM->G[i][j][k][0] = id;
            }         
        }        
    }
    else
    {
        cout << "Undefined boundary condition method." << endl;
        abort();
    }
    // Update neighour nodes list.
    UpdateLn(p0);
    cout << "finish" << endl;
}

inline void DELBM::ResetG()
{
    DomDEM->Lc.resize(0);
    // Reset flag for neighbor nodes.
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        for (size_t l=0; l<DomDEM->Lp[p]->Ln.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Ln[l][0];
            int j = DomDEM->Lp[p]->Ln[l][1];
            int k = DomDEM->Lp[p]->Ln[l][2];
            // Reset neighour nodes to fluid nodes (g=-1)
            if (i>=0 && i<=Nx && j>=0 && j<=Ny && k>=0 && k<=Nz)
            {
                if (DomLBM->G[i][j][k][0]<6. && DomLBM->G[i][j][k][0]>-1.)
                {
                    cout << "trying to rewrite G for boundary." << endl;
                    cout << i << " " << j << " " << k << endl;
                    abort();
                }
                DomLBM->G[i][j][k][0] = -1.;
            }
        }
    }
    Lr.resize(0);
}

inline void DELBM::UpdateGForIBB()
{
    ResetG();
    // vector < vector< vector< int > > >  Lc_threads(Nproc);
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        // Reset boundary nodes list.
        p0->Lb.clear();
        p0->Lq.clear();

        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Ln[l][0];
            int j = DomDEM->Lp[p]->Ln[l][1];
            int k = DomDEM->Lp[p]->Ln[l][2];

            Vector3d ind (i,j,k);
            double dis = (ind-p0->X).norm()-p0->R;

            bool mark = false;
            double id = p0->ID;

            if (dis<0.)     mark = true;
            // First check for possible crossing, sqrt(3) < 1.75 (for 3d)
            else if (dis<1.75)
            {
                VectorXd deltaQ(DomLBM->Q);
                for (int q=0; q<DomLBM->Q; ++q)
                {
                    deltaQ(q) = CalCrossPointQ(ind, DomLBM->E[DomLBM->Op[q]], p0->X, p0->R);
                    if (deltaQ(q)>-1.)      mark = true;
                }
                if (mark)
                {
                    #pragma omp critical
                    {
                        p0->Lb.push_back({i,j,k});
                        p0->Lq.push_back(deltaQ);
                    }
                    id += 0.5;
                    // Check distance between new boundary nodes and DEM particle surface at last time step to find refilling nodes.
                    if ((ind-p0->Xb).norm()-p0->R<0.)
                    {
                        #pragma omp critical
                        {
                            Lr.push_back({i,j,k});
                        }
                    }
                }
            }

            if (mark)
            {
                double g = DomLBM->G[i][j][k][0];
                if (g==-1.)     DomLBM->G[i][j][k][0] = id;
            }

            // if (mark)
            // {
            //     double g = DomLBM->G[i][j][k][0];
            //     if (g==-1.)     DomLBM->G[i][j][k][0] = id;
            //     else
            //     {
            //         if ((g-((int) g)) == 0.5 && (id-((int) id))==0)
            //         {
            //             DomLBM->G[i][j][k][0] = id;
            //         }

            //         int p = p0->ID;
            //         int q = (int) DomLBM->G[i][j][k][0];
            //         if (q==-2)
            //         {
            //             cout << "Contact" << endl;
            //             Vector3d n = ind-DomDEM->Lp[p]->X;
            //             // Overlapping distance
            //             double delta = DomDEM->Lp[p]->R+0.5-n.norm();
            //             if (delta>0)
            //             {
            //                 n.normalize();
            //                 Vector3d fc = n*pow(delta, 1.5);
            //                 DomDEM->Lp[p]->Fc += DomDEM->Lp[p]->Kn*fc;
            //             }
            //         }
            //         else
            //         {
            //             #pragma omp critical
            //             {
            //                 if (p!=q)   DomDEM->Lc.push_back({min(p,q), max(p,q)});
            //             }
            //         }
            //     }
            // }
        }
        // Update neighour nodes list.
        UpdateLn(p0);
    }
    sort( Lr.begin(), Lr.end() );
    Lr.erase( unique( Lr.begin(), Lr.end() ), Lr.end() );

    sort( DomDEM->Lc.begin(), DomDEM->Lc.end() );
    DomDEM->Lc.erase( unique( DomDEM->Lc.begin(), DomDEM->Lc.end() ), DomDEM->Lc.end() );

    // for (int i=0; i<DomDEM->Lc.size(); ++i)
    // {
    //     cout << "p= " << DomDEM->Lc[i][0] << endl;
    //     cout << "q= " << DomDEM->Lc[i][1] << endl;
    //     cout << "===============" << endl;
    // }
}

inline void DELBM::UpdateG()
{
    ResetG();

    // vector < vector< vector< int > > >  Lc_threads(Nproc);
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        // Reset boundary nodes list.
        p0->Lb.clear();
        p0->Ld.clear();
        p0->Li.clear();

        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Ln[l][0];
            int j = DomDEM->Lp[p]->Ln[l][1];
            int k = DomDEM->Lp[p]->Ln[l][2];

            Vector3d ind (i,j,k);
            double dis = (ind-p0->X).norm()-p0->R;

            bool mark = false;
            double id = p0->ID;

            if (dis<0.)     mark = true;
            // First check for possible crossing, sqrt(3) < 1.75 (for 3d)
            else if (dis<1.75)
            {
                for (int q=0; q<DomLBM->Q; ++q)
                {
                    double disn = (ind+DomLBM->E[q]-p0->X).norm()-p0->R;
                    if (disn<0.)
                    {
                        mark = true;
                        break;
                    }
                }
                if (mark)
                {
                    #pragma omp critical
                    {
                        Vector3d xi = Hi*(ind-p0->X).normalized()+ind;
                        p0->Lb.push_back({i,j,k});
                        p0->Ld.push_back(dis);
                        p0->Li.push_back(xi);
                    }
                    id += 0.5;
                    // Check distance between new boundary nodes and DEM particle surface at last time step to find refilling nodes.
                    if ((ind-p0->Xb).norm()-p0->R<0.)
                    {
                        #pragma omp critical
                        {
                            Lr.push_back({i,j,k});
                        }
                    }
                }
            }
            if (mark)
            {
                double g = DomLBM->G[i][j][k][0];
                if (g==-1.)     DomLBM->G[i][j][k][0] = id;
            }
        }
        // Update neighour nodes list.
        UpdateLn(p0);
    }
    sort( Lr.begin(), Lr.end() );
    Lr.erase( unique( Lr.begin(), Lr.end() ), Lr.end() );

    sort( DomDEM->Lc.begin(), DomDEM->Lc.end() );
    DomDEM->Lc.erase( unique( DomDEM->Lc.begin(), DomDEM->Lc.end() ), DomDEM->Lc.end() );
}

// inline void DELBM::UpdateG()
// {
//     DomDEM->Lc.resize(0);
//     // Reset flag for neighbor nodes.
//     // #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (size_t p=0; p<DomDEM->Lp.size(); ++p)
//     {
//         for (size_t l=0; l<DomDEM->Lp[p]->Ln.size(); ++l)
//         {
//             int i = DomDEM->Lp[p]->Ln[l][0];
//             int j = DomDEM->Lp[p]->Ln[l][1];
//             int k = DomDEM->Lp[p]->Ln[l][2];
//             // Reset neighour nodes to fluid nodes (g=-1)
//             if (i>=0 && i<=Nx && j>=0 && j<=Ny && k>=0 && k<=Nz)
//             {
//                 if (DomLBM->G[i][j][k][0]==-2.)
//                 {
//                     cout << "trying to rewrite G for boundary." << endl;
//                     cout << i << " " << j << " " << k << endl;
//                     abort();
//                 }
//                 DomLBM->G[i][j][k][0] = -1.;
//             }
//         }
//     }

//     Lr.resize(0);

//     // vector < vector< vector< int > > >  Lc_threads(Nproc);
//     // #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (size_t p=0; p<DomDEM->Lp.size(); ++p)
//     {
//         DEM_PARTICLE* p0 = DomDEM->Lp[p];
//         // Reset boundary nodes list.
//         p0->Lb.clear();
//         p0->Ld.clear();
//         p0->Li.clear();

//         for (size_t l=0; l<p0->Ln.size(); ++l)
//         {
//             int i = DomDEM->Lp[p]->Ln[l][0];
//             int j = DomDEM->Lp[p]->Ln[l][1];
//             int k = DomDEM->Lp[p]->Ln[l][2];

//             Vector3d ind (i,j,k);
//             double dis = (ind-p0->X).norm()-p0->R;

//             bool mark = false;
//             double id = p0->ID;

//             if (dis<0.)     mark = true;
//             // First check for possible crossing, sqrt(3) < 1.75 (for 3d)
//             else if (dis<1.75)
//             {
//                 for (int q=0; q<DomLBM->Q; ++q)
//                 {
//                     double disn = (ind+DomLBM->E[q]-p0->X).norm()-p0->R;
//                     if (disn<0.)
//                     {
//                         mark = true;
//                         break;
//                     }
//                 }
//                 if (mark)
//                 {
//                     #pragma omp critical
//                     {
//                         Vector3d xi = Hi*(ind-p0->X).normalized()+ind;
//                         p0->Lb.push_back({i,j,k});
//                         p0->Ld.push_back(dis);
//                         p0->Li.push_back(xi);
//                     }
//                     id += 0.5;
//                 	// Check distance between new boundary nodes and DEM particle surface at last time step to find refilling nodes.
//                 	if ((ind-p0->Xb).norm()-p0->R<0.)
//                     {
//                         #pragma omp critical
//                         {
//                             Lr.push_back({i,j,k});
//                         }
//                     }
//                 }
//             }
//             if (mark)
//             {
//                 double g = DomLBM->G[i][j][k][0];
//                 if (g==-1.)     DomLBM->G[i][j][k][0] = id;
//             }
//             // if (mark)
//             // {
//             //     double g = DomLBM->G[i][j][k][0];
//             //     if (g==-1.)     DomLBM->G[i][j][k][0] = id;
//             //     else
//             //     {
//             //         if ((g-((int) g)) == 0.5 && (id-((int) id))==0)
//             //         {
//             //             DomLBM->G[i][j][k][0] = id;
//             //         }

//             //         int p = p0->ID;
//             //         int q = (int) DomLBM->G[i][j][k][0];
//             //         if (q==-2)
//             //         {
//             //             cout << "Contact" << endl;
//             //             Vector3d n = ind-DomDEM->Lp[p]->X;
//             //             // Overlapping distance
//             //             double delta = DomDEM->Lp[p]->R+0.5-n.norm();
//             //             if (delta>0)
//             //             {
//             //                 n.normalize();
//             //                 Vector3d fc = n*pow(delta, 1.5);
//             //                 DomDEM->Lp[p]->Fc += DomDEM->Lp[p]->Kn*fc;
//             //             }
//             //         }
//             //         else
//             //         {
//             //             #pragma omp critical
//             //             {
//             //                 if (p!=q)   DomDEM->Lc.push_back({min(p,q), max(p,q)});
//             //             }
//             //         }
//             //     }
//             // }
//         }
//         // Update neighour nodes list.
//         UpdateLn(p0);
//     }
//     sort( Lr.begin(), Lr.end() );
//     Lr.erase( unique( Lr.begin(), Lr.end() ), Lr.end() );

//     sort( DomDEM->Lc.begin(), DomDEM->Lc.end() );
//     DomDEM->Lc.erase( unique( DomDEM->Lc.begin(), DomDEM->Lc.end() ), DomDEM->Lc.end() );

//     // for (int i=0; i<DomDEM->Lc.size(); ++i)
//     // {
//     //     cout << "p= " << DomDEM->Lc[i][0] << endl;
//     //     cout << "q= " << DomDEM->Lc[i][1] << endl;
//     //     cout << "===============" << endl;
//     // }
// }

inline void DELBM::UpdateGForVAM()
{
    // Reset flag for neighbor nodes.
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];
            // If the node is not domain boundary nodes
            if (DomLBM->G[i][j][k][0] > 5.)
            {
                DomLBM->G[i][j][k][0] = -1.;
                DomLBM->G[i][j][k][1] =  0.;
            }
            DomLBM->G[i][j][k][2] = -1.;
            DomLBM->G[i][j][k][3] =  0.;
        }      
    }
    // vector < vector< vector< int > > >  Lc_threads(Nproc);
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        // Reset boundary nodes list.
        p0->Lb.clear();
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];
            double tempG = CalGdistance(i,j,k,p0);
            // double tempG = 0;
            // if (Gmethod==0)            tempG = CalGdistance(i,j,k,p0);
            // else if (Gmethod==1)       tempG = CalGvolumeSphere(i,j,k,p0, 20);
            // else
            // {
            //     cout << "wrong Gmethod!" << endl;
            //     abort();
            // }
            if (tempG>0)
            {
                if (tempG<1.)
                {
                    p0->Lb.push_back({i,j,k});
                }
                if (DomLBM->G[i][j][k][0]==-1.)
                {
                    DomLBM->G[i][j][k][0] = p;
                    DomLBM->G[i][j][k][1] = tempG;
                }
                else if (DomLBM->G[i][j][k][2]==-1.)
                {
                    DomLBM->G[i][j][k][2] = p;
                    DomLBM->G[i][j][k][3] = tempG;

                    // int q = (int) DomLBM->G[i][j][k][0];
                    // Lc_threads[omp_get_thread_num()].push_back({min(p,q), max(p,q)});
                }
                else
                {
                    // int q = (int) DomLBM->G[i][j][k][2];
                    // Lc_threads[omp_get_thread_num()].push_back({min(p,q), max(p,q)});  
                }
            }
        }
        // Update neighour nodes list.
        UpdateLn(p0);
    }
}

inline void DELBM::UpdateLn(DEM_PARTICLE* p0)
{
    p0->Ln.clear();
    for (size_t l=0; l<p0->Lb.size(); ++l)
    {
        for (size_t q=0; q<DomLBM->Ne.size(); ++q)
        {
            int i = p0->Lb[l][0]+DomLBM->Ne[q][0];
            int j = p0->Lb[l][1]+DomLBM->Ne[q][1];
            int k = p0->Lb[l][2]+DomLBM->Ne[q][2];

            p0->Ln.push_back({i,j,k});
            // int i = p0->Lb[l][0]+DomLBM->Ne[q][0];
            // int j = p0->Lb[l][1]+DomLBM->Ne[q][1];
            // int k = p0->Lb[l][2]+DomLBM->Ne[q][2];

            // if (Periodic[0])    i = (i+Nx+1)%(Nx+1);
            // if (Periodic[1])    j = (j+Ny+1)%(Ny+1);
            // if (Periodic[2])    k = (k+Nz+1)%(Nz+1);

            // if (i>=0 && i<=Nx && j>=0 && j<=Ny && k>=0 && k<=Nz)
            // {
            //     p0->Ln.push_back({i,j,k});
            //     if (DomLBM->G[i][j][k][0]!=-2.)
            //     {
            //         #pragma omp critical
            //         {
            //             p0->Ln.push_back({i,j,k});
            //         }
            //     }
            // }
            // cout << p0->Lb[l][0] << " " << p0->Lb[l][1] << " " << p0->Lb[l][2] << " " << endl;
            // cout << i << " " << j << " " << k << " " << endl;
            // cout << "(p0->Lb[l][0]+DomLBM->Ne[q][0]+PeriodicIndex[0]*(Nx+1))= " << (p0->Lb[l][0]+DomLBM->Ne[q][0]+PeriodicIndex[0]*(Nx+1)) << endl;
            // cout << "(PeriodicIndex[0]*Nx+1)= " << (PeriodicIndex[0]*Nx+1) << endl;
            // abort();
        }
    }
    sort( p0->Ln.begin(), p0->Ln.end() );
    p0->Ln.erase( unique( p0->Ln.begin(), p0->Ln.end() ), p0->Ln.end() );
}

inline double DELBM::CalGdistance(int i, int j, int k, DEM_PARTICLE* p0)
{
    Vector3d xc (i,j,k);
    double dis = (xc-p0->X).norm()-p0->R;
    double g = 0.;

    if (dis<-0.5)  g = 1.;
    else if (dis<0.5)   g = 0.5-dis;
    return g;
}

inline double DELBM::CalGvolumeSphere(int i, int j, int k, DEM_PARTICLE* p0, int n)
{
    int nc = 0.;
    double v = 0.;

    Vector3d xl (i,j,k);
    double dis = (xl-p0->X).norm()-p0->R;

    if (dis<-0.866)     v = 1.;
    else if (dis>0.866) v = 0.;
    else
    {
        Vector3d xc (0.,0.,0.);

        for (int ni=0; ni<n; ++ni)
        for (int nj=0; nj<n; ++nj)
        for (int nk=0; nk<n; ++nk)
        {
            xc(0) = i-0.5+(ni+0.5)/n;
            xc(1) = j-0.5+(nj+0.5)/n;
            xc(2) = k-0.5+(nk+0.5)/n;

            if ((xc-p0->X).norm()<p0->R)    nc++;
        }

        v = ((double) nc)/((double) n*n*n);
    }
    return v;
}

// If not crossing, q is set to -1.
inline double DELBM::CalCrossPointQ(Vector3d& a, Vector3d& v, Vector3d& c, double r)
{
    double q = -1.;

    Vector3d ac = a-c;
    double A = v.squaredNorm();
    double B = 2.*v.dot(ac);
    double C = ac.squaredNorm()-r*r;

    double D = B*B-4.*A*C;

    if (D>=0.)
    {
        D = sqrt(D);

        double q0 = 0.5*(-B+D)/A;
        double q1 = 0.5*(-B-D)/A;

        if      (q0>=0 && q0<=1)    q = q0;
        else if (q1>=0 && q1<=1)    q = q1;
    }

    return q;
}

inline void DELBM::Refill()
{
    for (size_t l=0; l<Lr.size(); ++l)
    {
            int i = Lr[l][0];
            int j = Lr[l][1];
            int k = Lr[l][2];

            Vector3d ind (i,j,k);

            int p = (int) DomLBM->G[i][j][k][0];

            DEM_PARTICLE* p0 = DomDEM->Lp[p];

            VectorXd feqi(DomLBM->Q);

            double rhoi = 0.;
            int addRho = 0;

            // Vector3d vi = p0->V + p0->W.cross(ind - p0->X);

            for (int q=1; q<DomLBM->Q; ++q)
            {
                int in = i + ((int) DomLBM->E[q](0));
                int jn = j + ((int) DomLBM->E[q](1));
                int kn = k + ((int) DomLBM->E[q](2));

                double g = DomLBM->G[in][jn][kn][0];
                // double gt = DomLBM->Gt[in][jn][kn];

                // Not from fluid nodes
                if (g==-1. || (g-((int) g))==0.5)
                {

                    Vector3d aa = -DomLBM->E[q] + p0->X - p0->Xb;
                    Vector3d bb = ind - p0->X;

                    double A = aa.squaredNorm();
                    double B = 2.*aa.dot(bb);
                    double C = bb.squaredNorm() - p0->R*p0->R;
                    double D = B*B-4.*A*C;

                    double delta = -1.;

                    if (D>=0.)
                    {
                        D = sqrt(D);

                        double delta0 = 0.5*(-B+D)/A;
                        double delta1 = 0.5*(-B-D)/A;

                        if      (delta0>=0 && delta0<=1)    delta = delta0;
                        else if (delta1>=0 && delta1<=1)    delta = delta1;
                    }

                    Vector3d xk = p0->Xb + (1.-delta)*(p0->X - p0->Xb);

                    Vector3d n = ind-delta*DomLBM->E[q]-xk;
                    Vector3d vp = p0->V + p0->W.cross(n);

                    DomLBM->F[i][j][k](q) =  DomLBM->Ft[in][jn][kn](DomLBM->Op[q]) + 6.*DomLBM->W[q]*DomLBM->Rho[i][j][k]*DomLBM->E[q].dot(vp);
                }
                else
                {
                    rhoi += DomLBM->Rho[in][jn][kn];
                    addRho++;
                }
            }

            // rhoi /= addRho;

            // DomLBM->CalFeqC(feqi, rhoi, vi);
            // DomLBM->F[i][j][k](0) = feqi(0);
            // for (int i=0; i<DomLBM->Ne.size(); ++i)
            // {
            //     int ine = i+DomLBM->Ne[q][0];
            //     int jne = j+DomLBM->Ne[q][1];
            //     int kne = k+DomLBM->Ne[q][2];

            //     if (DomLBM->G[ine][jne][kne]!=-1 && DomLBM->G[ine][jne][kne]!=p)
            //     {

            //     }                
            // }
    }
}

inline void DELBM::NormRefill()
{
    for (size_t l=0; l<Lr.size(); ++l)
    {
        int i = Lr[l][0];
        int j = Lr[l][1];
        int k = Lr[l][2];

        Vector3d ind (i,j,k);

        int p = (int) DomLBM->G[i][j][k][0];

        DEM_PARTICLE* p0 = DomDEM->Lp[p];

        Vector3d n = (ind-p0->X).normalized();

        double alpha = -1.0e12;
        int ord = 1;

        for (int q=1; q<DomLBM->Q; ++q)
        {
            double dot = n.dot(DomLBM->E[q])/DomLBM->E[q].norm();

            if (dot>alpha)
            {
                alpha = dot;
                ord++;
            }
        }

        int in = i + ((int) DomLBM->E[ord](0));
        int jn = j + ((int) DomLBM->E[ord](1));
        int kn = k + ((int) DomLBM->E[ord](2));

        int inn = i + 2*((int) DomLBM->E[ord](0));
        int jnn = j + 2*((int) DomLBM->E[ord](1));
        int knn = k + 2*((int) DomLBM->E[ord](2));

        int innn = i + 3*((int) DomLBM->E[ord](0));
        int jnnn = j + 3*((int) DomLBM->E[ord](1));
        int knnn = k + 3*((int) DomLBM->E[ord](2));           

        // DomLBM->F[i][j][k] = 3.*DomLBM->F[in][jn][kn]-3.*DomLBM->F[inn][jnn][knn]+DomLBM->F[innn][jnnn][knnn];

        DomLBM->F[i][j][k](0) = 3.*DomLBM->F[in][jn][kn](0)-3.*DomLBM->F[inn][jnn][knn](0)+DomLBM->F[innn][jnnn][knnn](0);

        for (int q=1; q<DomLBM->Q; ++q)
        {
            int ip = (i - (int) DomLBM->E[q][0]);
            int jp = (j - (int) DomLBM->E[q][1]);
            int kp = (k - (int) DomLBM->E[q][2]);

            if (DomLBM->G[ip][jp][kp][0]==0.)
            {
                DomLBM->F[i][j][k](q) = 3.*DomLBM->F[in][jn][kn](q)-3.*DomLBM->F[inn][jnn][knn](q)+DomLBM->F[innn][jnnn][knnn](q);
            }
        }
    }
}

inline void DELBM::EqRefill()
{
    for (size_t l=0; l<Lr.size(); ++l)
    {
        int i = Lr[l][0];
        int j = Lr[l][1];
        int k = Lr[l][2];
        // cout << "i= " << i << " j= " << j << " k= " << k << endl;
        Vector3d ind (i,j,k);
        int p = (int) DomLBM->G[i][j][k][0];
        // cout << "p= " << p << endl;
        // cout << "DomLBM->G[i][j][k][0]=  " << DomLBM->G[i][j][k][0] << endl;
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        Vector3d n = ind-p0->X;
        Vector3d vp = p0->V + p0->W.cross(n);
        VectorXd feq (DomLBM->Q);
        DomLBM->CalFeqC(feq, DomLBM->Rho0, vp);
        DomLBM->F[i][j][k](0) = feq(0);
        for (int q=1; q<DomLBM->Q; ++q)
        {
            int ip = (i - (int) DomLBM->E[q][0]);
            int jp = (j - (int) DomLBM->E[q][1]);
            int kp = (k - (int) DomLBM->E[q][2]);

            if (DomLBM->G[ip][jp][kp][0]==6.)
            {
                DomLBM->F[i][j][k](q) = feq(q);
            }
        }
    }
}
// Apply Volume Average Method
// inline void DELBM::ApplyVAM()
// {
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (int n=0; n<DomLBM->Ncell; ++n)
//     {
//         int i, j, k;
//         DomLBM->FindIndex(n, i, j, k);
//         DomLBM->CollideSRTLocal(i, j, k);
//         if (DomLBM->G[i][j][k][0]>5.)
//         {
//             double g = DomLBM->G[i][j][k][1];
//             // g = g*(DomLBM->Tau-0.5)/(1.-g + DomLBM->Tau-0.5);
//             size_t p = DomLBM->G[i][j][k][0];
//             DEM_PARTICLE* p0 = DomDEM->Lp[p];
//             double rhos = p0->Rho;
//             Vector3d ind (i,j,k);
//             Vector3d n = (ind-p0->X).normalized()*p0->R;
//             if (g==2.)
//             {

//             }

//             Vector3d vp = p0->V + p0->W.cross(n);
//             Vector3d fh (0.,0.,0.);
//             DomLBM->VAM(i,j,k,g,rhos,vp,fh);
//             Vector3d th = n.cross(fh);

//             p0->Fh += fh;
//             p0->Th += th;
//         }
//     }
// }

// Apply Volume Average Method
inline void DELBM::ApplyPSM()
{
    DomDEM->CMap.clear();
    // Reset flag for neighbor nodes.
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];

            int ic = (i+(Nx+1))%(Nx+1);
            int jc = (j+(Ny+1))%(Ny+1);
            int kc = (k+(Nz+1))%(Nz+1);
            
            // If the node is not domain boundary nodes
            if (DomLBM->G[ic][jc][kc][0] > 5.)
            {
                DomLBM->G[ic][jc][kc][0] = -1.;
                DomLBM->G[ic][jc][kc][1] =  0.;
            }
            DomLBM->G[i][j][k][2] = -1.;
            DomLBM->G[i][j][k][3] =  0.;
        }
    }

    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        p0->Lb.clear();
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];

            int ic = (i+(Nx+1))%(Nx+1);
            int jc = (j+(Ny+1))%(Ny+1);
            int kc = (k+(Nz+1))%(Nz+1);

            // double g = CalGdistance(i,j,k,p0);
            double g = CalGvolumeSphere(i,j,k,p0, 20);
            bool sameCell = false;
            // bool contact = false;
            if (g>0. && g<1.)
            {
                if (DomLBM->G[ic][jc][kc][0] == -1.)  
                {
                    DomLBM->G[ic][jc][kc][0] = p;
                    DomLBM->G[ic][jc][kc][1] = g;
                }
                else
                {
                    size_t q = DomLBM->G[ic][jc][kc][0];
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
                        if (!DomDEM->CMap[key])
                        {
                            DomDEM->Lc.push_back({min0, max0});
                            DomDEM->CMap[key] = true;
                        }
                        sameCell = true;
                    }
                }
                p0->Lb.push_back({i,j,k});
                Vector3d ind (i,j,k);
                Vector3d n = (ind-p0->X).normalized()*p0->R;
                Vector3d vp = p0->V + p0->W.cross(n);
                Vector3d fh (0.,0.,0.);
                // DomLBM->VAM(ic,jc,kc,g,rhos,vp,fh);
                // Vector3d th = n.cross(fh);
                if (!sameCell)
                {
                    DomLBM->PSM(ic,jc,kc,g,vp,fh);
                    Vector3d th = n.cross(fh); 
                    p0->Fh += fh;
                    p0->Th += th;                   
                }
                // p0->Fh += fh;
                // p0->Th += th;                
            }
        }
        UpdateLn(p0);      
    }
}

// Apply Volume Average Method
inline void DELBM::ApplyVAM()
{
    DomDEM->CMap.clear();
    // Reset flag for neighbor nodes.
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];

            int ic = i;
            int jc = j;
            int kc = k;

            bool out = false;
            if (Periodic[0])            ic = (i+(Nx+1))%(Nx+1);
            else if (ic<0 || ic>Nx)     out=true;
            if (Periodic[1])            jc = (j+(Ny+1))%(Ny+1);
            else if (jc<0 || jc>Ny)     out=true;
            if (Periodic[2])            kc = (k+(Nz+1))%(Nz+1);
            else if (kc<0 || kc>Nz)     out=true;
            
            if (!out)
            {
                if (DomLBM->Flag[ic][jc][kc].size()>0)
                {
                    vector<size_t> flagtmp;
                    flagtmp.resize(0);
                    for (size_t m=0; m<DomLBM->Flag[ic][jc][kc].size(); ++m)
                    {
                        size_t ind = DomLBM->Flag[ic][jc][kc][m];
                        if (ind<6)
                        {
                            flagtmp.push_back(ind);
                        }
                    }
                    DomLBM->Flag[ic][jc][kc] = flagtmp;
                }
            }
        }
    }
    // cout << "reset finish" << endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        // cout << "p= " << p << endl;
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        p0->Lb.clear();
        // cout << "p0->Lb.clear();" << endl;
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];

            int ic = i;
            int jc = j;
            int kc = k;

            bool out = false;
            if (Periodic[0])            ic = (i+(Nx+1))%(Nx+1);
            else if (ic<0 || ic>Nx)     out=true;
            if (Periodic[1])            jc = (j+(Ny+1))%(Ny+1);
            else if (jc<0 || jc>Ny)     out=true;
            if (Periodic[2])            kc = (k+(Nz+1))%(Nz+1);
            else if (kc<0 || kc>Nz)     out=true;

            if (!out)
            {
                double g = CalGdistance(i,j,k,p0);
                Vector3d cell (i,j,k);
                double dis = (cell-p0->X).norm()-p0->R;
                // cout << "CalGdistance(i,j,k,p0);" << endl;
                // double g = CalGvolumeSphere(i,j,k,p0, 20);
                // double bn = g*(DomLBM->Tau-0.5)/(1.-g + DomLBM->Tau-0.5);
                double bn = g;
                if (g>0. && g<1.)
                // if (dis<0.87 && g>0.)
                {
                    #pragma omp critical
                    {
                        if (DomLBM->Flag[ic][jc][kc].size()>0)
                        {
                            for (size_t m=0; m<DomLBM->Flag[ic][jc][kc].size(); ++m)
                            {
                                size_t q = DomLBM->Flag[ic][jc][kc][m];
                                if (q<6 && DomDEM->Periodic[q/2])
                                {
                                    p0->crossing[q/2] = true;
                                    p0->crossingFlag = true;
                                }
                                else
                                {
                                    size_t min0 = min(p,q);
                                    size_t max0 = max(p,q);
                                    size_t key = Key(min0,max0);
                                    if (!DomDEM->CMap[key])
                                    {
                                        DomDEM->Lc.push_back({min0, max0});
                                        DomDEM->CMap[key] = true;
                                    }
                                }       
                            }
                        }
                        // DomLBM->Flag[ic][jc][kc].push_back(p);                    
                    }

                    p0->Lb.push_back({i,j,k});
                    // cout << "p0->Lb.push_back({i,j,k});" << endl;
                    double rhos = p0->Rho;
                    Vector3d ind (i,j,k);
                    Vector3d n = (ind-p0->X).normalized()*p0->R;
                    Vector3d vp = p0->V + p0->W.cross(n);
                    Vector3d fh (0.,0.,0.);
                    DomLBM->VAM(ic,jc,kc,bn,rhos,vp,fh);
                    // if (p0->ID==1000)
                    // {
                    //     cout << "bn= " << bn << endl;
                    //     cout << "fh= " << fh.transpose() << endl;
                    //     cout << "vp= " << vp.transpose() << endl;
                    // }
                    // cout << "DomLBM->VAM(ic,jc,kc,bn,rhos,vp,fh);" << endl;
                    Vector3d th = n.cross(fh);
                    p0->Fh += fh;
                    p0->Th += th;   
                    // cout << "p0->Fh += fh;" << endl;             
                }
                if (dis<0.87)
                {
                    #pragma omp critical
                    {
                        DomLBM->Flag[ic][jc][kc].push_back(p); 
                    }
                }
            }
        }
        UpdateLn(p0);
        // cout << "UpdateLn(p0);" << endl;  
    }
}

inline void DELBM::UpdateBoxGlobal()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=6; i<DomDEM->Lp.size(); ++i)
    {
        DomDEM->Lp[i]->UpdateBox(D);
    }    
    for (size_t i=0; i<DomDEM->Lg.size(); ++i)
    {
        DEM_PARTICLE* g0 = DomDEM->Lg[i];
        g0->UpdateBox(D);
        for (size_t l=0; l<g0->Lp.size(); ++l)
        {
            size_t p = g0->Lp[l];
            DomDEM->Lp[p]->UpdateBox(D);
        }
    }
}

inline void DELBM::MoveRW()
{
    #pragma omp parallel for schedule(static) num_threads(DomRWM->Nproc)
    for (size_t p=0; p<DomRWM->Lp.size(); ++p)
    {
        RWM_PARTICLE* p0 = DomRWM->Lp[p];
        // Store the position before move
        p0->Xb = p0->X;
        // convection velocity
        Vector3d vc = DomLBM->InterpolateV(p0->X);

        for(int d=0; d<D; ++d)
        {
            p0->X(d) += vc(d)*DomRWM->Dt + DomRWM->GetNormalD()*sqrt(2*DomRWM->Dc*DomRWM->Dt);
            // if (Lp[p]->X(d)<0 || Lp[p]->X(d)>DomSize[d]) Lrm.push_back(p);
        }
        Vector3d xm = p0->X;
        p0->Vc = p0->X;
        // find which cell this rw particles belongs
        int i = round(p0->X(0));
        int j = round(p0->X(1));
        int k = round(p0->X(2));
        i = (i+Nx+1)%(Nx+1);
        j = (j+Ny+1)%(Ny+1);
        k = (k+Nz+1)%(Nz+1);
        // if its ocuppied by dem
        // bool did = false;
        if (DomLBM->Flag[i][j][k].size()>0)
        {
            // list of real particles
            vector<size_t> li;
            vector<Vector3d> lx;
            // the primary particle used to reflect rw
            size_t ind0 = 0;
            Vector3d xp0;
            // shortest distane to determine primary particle
            double shortestDis = 1.0e12;
            for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
            {
                size_t ind = DomLBM->Flag[i][j][k][m];
                DEM_PARTICLE* demP = DomDEM->Lp[ind];
                // for real particles
                if (ind>5)
                {
                    double dis = 1.0e12;
                    Vector3d xpt;
                    // if dem is crossing BC then check the shortest distance for every mirror positions
                    if (demP->crossingFlag)
                    {
                        for (size_t i=0; i<demP->Xmir.size(); ++i)
                        {
                            double dist = (demP->Xmir[i] - p0->X).norm()-demP->R;
                            if (dist<dis)
                            {
                                xpt = demP->Xmir[i];
                                dis = dist;
                            }
                        }
                    }
                    else
                    {
                        xpt = demP->X;
                        dis = (demP->X-p0->X).norm()-demP->R;
                    }
                    li.push_back(ind);
                    lx.push_back(xpt);

                    if (dis<shortestDis)
                    {
                        ind0 = ind;
                        shortestDis = dis;
                        xp0 = xpt;
                    }
                }
            }
            // check for collision with walls
            for (size_t m=0; m<3; ++m)
            {
                size_t ind = DomLBM->Flag[i][j][k][m];
                if (ind<6 && !Periodic[ind/2])
                {
                    li.push_back(ind);
                    double dis = (p0->X(ind/2)-DomDEM->Lp[ind]->X(ind/2))*DomDEM->Lp[ind]->Normal(ind/2);
                    if (dis<shortestDis)
                    {
                        ind0 = ind;
                        shortestDis = dis;
                        xp0 = DomDEM->Lp[ind]->X;
                    }
                }
            }
            // if do contact with dem
            if (shortestDis<0.)
            {
                DEM_PARTICLE* demP0 = DomDEM->Lp[ind0];
                Vector3d xc, n;
                double r = demP0->R;
                // find contact point and normal for boundaries
                if (ind0<6)
                {
                    size_t dirc = ind0/2;
                    n = demP0->Normal;
                    // n.setZero();
                    // n(dirc) = pow(-1, ind0%2);
                    double xbd = p0->Xb(dirc)-demP0->X(dirc);
                    double xd = p0->X(dirc)-demP0->X(dirc);
                    double k = abs(xbd/(xd-xbd));
                    xc = p0->Xb + k*(p0->X-p0->Xb);
                }
                // find contact point and normal for real particles
                else
                {
                    Vector3d xp0b = xp0+(demP0->Xbr-demP0->X);
                    DomRWM->ContactPointWithSphere(p0, xp0, xp0b, r, xc, n);
                }
                        
                DomRWM->Reflection(p0, n, xc);
                // if still within dem particle then move to surface along normal direction
                if ((p0->X-demP0->X).squaredNorm()<r*r)
                {
                    p0->X = xp0+(r+1.0e-12)*(p0->X-xp0).normalized();
                }
                // double check if rm particle is within other particles after reflection
                if (li.size()>1)
                {
                    bool inside = false;
                    for (size_t m=0; m<li.size(); ++m)
                    {
                        size_t ind = li[m];
                        double dis;
                        if (ind>5)  dis = (lx[m]-p0->X).norm()-DomDEM->Lp[ind]->R;
                        else        dis = (p0->X(ind/2)-DomDEM->Lp[ind]->X(ind/2))*DomDEM->Lp[ind]->Normal(ind/2);
                        if (dis<0.)
                        {
                            inside = true;
                            break;
                        }
                    }
                    if (inside)
                    {
                        // the position of rw after convection without diffusion
                        Vector3d xa = p0->Xb + vc*DomRWM->Dt;
                        // cout << "xa: " << xa.transpose() << endl;
                        for (size_t m=0; m<li.size(); ++m)
                        {
                            for (size_t n=0; n<1.0e12; ++n)
                            {
                                Vector3d xcor = xa+sqrt(2*DomRWM->Dc*DomRWM->Dt)*DomRWM->GetRandomPointinUnitSphere();
                                // cout << "xa: " << xa.transpose() << endl;
                                bool ok = true;
                                for (size_t m=0; m<li.size(); ++m)
                                {
                                    size_t ind = li[m];
                                    double dis;
                                    if (ind>5) dis = (lx[m]-xcor).norm()-DomDEM->Lp[ind]->R;
                                    else dis = (xcor(ind/2)-DomDEM->Lp[ind]->X(ind/2))*pow(-1, ind%2);
                                    if (dis<0.)
                                    {
                                        ok = false;
                                        break;
                                    }
                                    // cout << "lx[m]: " << lx[m].transpose() << endl;
                                    // cout << "dis: " << (xcor-lx[m]).norm() << endl;
                                }
                                if (ok)
                                {
                                    // cout << "xa is outside" << endl;
                                    p0->X = xcor;
                                    // cout << "xcor: " << xcor.transpose() << endl;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // // if (li)
            // for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
            // // for (size_t m=0; m<1; ++m)
            // {
            //     size_t ind = DomLBM->Flag[i][j][k][m];
            //     if (ind>=DomDEM->Lp.size())
            //     {
            //         cout << "wrong ind: " << ind << endl;
            //         cout << "DomLBM->Flag[i][j][k][m]: " << DomLBM->Flag[i][j][k][m] << endl;
            //         cout << DomDEM->Lp.size() << endl;
            //         cout << i << " " << j << " " << k << endl;
            //     }
            //     if (ind>=6)
            //     {
            //         double dis2 = (DomDEM->Lp[ind]->X-p0->X).squaredNorm();
            //         double r = DomDEM->Lp[ind]->R;
            //         double r2 = DomDEM->Lp[ind]->R*DomDEM->Lp[ind]->R;
            //         if (dis2<r2)
            //         {
            //             did = true;
            //             // double disb2 = (DomDEM->Lp[ind]->Xbr-p0->Xb).squaredNorm();
            //             // if (disb2<r2)
            //             // {
            //             //     cout << " its inside already" << endl;
            //             //     cout << sqrt(disb2) << endl;
            //             //     cout << sqrt(dis2) << endl;
            //             //     cout << "p:" << p << endl;
            //             //     cout << "DomDEM->Lp[ind]->Xbr: " << DomDEM->Lp[ind]->Xbr.transpose() << endl;
            //             //     cout << "DomDEM->Lp[ind]->X: " << DomDEM->Lp[ind]->X.transpose() << endl;
            //             //     cout << "p0->Xb: " << p0->Xb.transpose() << endl;
            //             //     cout << "p0->X: " << p0->X.transpose() << endl;
            //             //     abort();
            //             // }
            //             Vector3d xc, n;
            //             DomRWM->ContactPointWithSphere(p0, DomDEM->Lp[ind]->X, DomDEM->Lp[ind]->Xbr, DomDEM->Lp[ind]->R, xc, n);                        
            //             DomRWM->Reflection(p0, n, xc);
            //             if ((p0->X-DomDEM->Lp[ind]->X).squaredNorm()<r2)
            //             {
            //                 p0->X = DomDEM->Lp[ind]->X+(r+1.0e-12)*(p0->X-DomDEM->Lp[ind]->X).normalized();
            //             }
            //             // p0->Xb = xc;
            //             // if ((p0->X-DomDEM->Lp[ind]->X).squaredNorm()<r2)
            //             // {
            //             //     cout << "inside.........." << endl;
            //             //     abort();
            //             // }
            //         }
            //     }
            // }
            // bool stop = false;
            // for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
            // {
            //     size_t ind = DomLBM->Flag[i][j][k][m];
            //     if (ind>=6)
            //     {
            //         double dis2 = (DomDEM->Lp[ind]->X-p0->X).squaredNorm();
            //         double r = DomDEM->Lp[ind]->R;
            //         double r2 = DomDEM->Lp[ind]->R*DomDEM->Lp[ind]->R;
            //         if (dis2<r2)
            //         {
            //             cout << "still inside" << endl;
            //             // cout << "p0->Xb: " << p0->Xb.transpose() << endl;
            //             // cout << "p0->Vc: " << p0->Vc.transpose() << endl;
            //             stop = true;
            //             cout << "inside of ind: " << ind << endl;
            //         }
            //     }
            // }
            // if (stop)
            // {
            //     vector<size_t> li;
            //     vector<Vector3d> lx;
            //     vector<double> lr;
            //     for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
            //     {
            //         size_t ind = DomLBM->Flag[i][j][k][m];
            //         if (ind>5)
            //         {
            //             li.push_back(ind);
            //             lx.push_back(DomDEM->Lp[ind]->X);
            //             lr.push_back(DomDEM->Lp[ind]->R);
            //         }
            //         cout << "ind: " << ind << endl;
            //     }
            //     // abort();
            //     Vector3d xa = p0->Xb + vc*DomRWM->Dt;

            //     for (size_t n=0; n<1.0e12; ++n)
            //     {
            //         Vector3d xcor = xa+sqrt(2*DomRWM->Dc*DomRWM->Dt)*DomRWM->GetRandomPointinUnitSphere();
            //         cout << "xa: " << xa.transpose() << endl;
            //         bool ok = true;
            //         for (size_t m=0; m<li.size(); ++m)
            //         {
            //             if ((xcor-lx[m]).norm()<lr[m])  ok = false;
            //             cout << "lx[m]: " << lx[m].transpose() << endl;
            //             cout << "dis: " << (xcor-lx[m]).norm() << endl;
            //         }
            //         if (ok)
            //         {
            //             // cout << "xa is outside" << endl;
            //             p0->X = xcor;
            //             cout << "xcor: " << xcor.transpose() << endl;
            //             // abort();
            //             break;
            //         }
            //         // abort();
            //     }
            // }
            size_t bc = DomLBM->Flag[i][j][k][0];
            if (bc<6)
            {
                if      (p0->X(0)>Nx)   p0->X(0) = p0->X(0)-Nx-1;
                else if (p0->X(0)<0.)   p0->X(0) = p0->X(0)+Nx+1;
                if      (p0->X(1)>Ny)   p0->X(1) = p0->X(1)-Ny-1;
                else if (p0->X(1)<0.)   p0->X(1) = p0->X(1)+Ny+1;
                if      (p0->X(2)>Nz)   p0->X(2) = p0->X(2)-Nz-1;
                else if (p0->X(2)<0.)   p0->X(2) = p0->X(2)+Nz+1;
            }
        }
        // double dis = (p0->X-DomDEM->Lp[6]->X).norm()-DomDEM->Lp[6]->R;
        // if (dis<0)
        // {
        //     if (did)
        //     {
        //         cout << "did but still insdie" << endl;
        //     }
        //     else
        //     {
        //         cout << "dont did and insdide" << endl;
        //     }
        //     // abort();
        // }
        // did = true;
        // if (did)
        // {
        //     // double dis = (p0->X-DomDEM->Lp[6]->X).norm()-DomDEM->Lp[6]->R;
        //     if (dis<0)
        //     {
        //         cout << "shit happened" << endl;
        //         cout << DomLBM->Flag[i][j][k].size() << endl;
        //         Vector3d cell (i,j,k);
        //         cout << "cell: " << cell.transpose() << endl;
        //         cout << "p0->X: " << p0->X.transpose() << endl;
        //         cout << (cell-DomDEM->Lp[6]->X).norm()-7.5 << endl;
        //         cout << (cell-DomDEM->Lp[6]->Xbr).norm()-7.5 << endl;  
        //         for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
        //         {
        //             size_t ind = DomLBM->Flag[i][j][k][m];
        //             if (ind>=6)
        //             {
        //                 double dis2 = (DomDEM->Lp[ind]->X-p0->X).squaredNorm();
        //                 double r2 = DomDEM->Lp[ind]->R*DomDEM->Lp[ind]->R;
        //                 cout << dis2-r2 << endl;
        //             }
        //         }
        //         abort();
        //     }
        // }
        if (p0->X(0)>Nx || p0->X(0)<0 || p0->X(1)>Ny || p0->X(1)<0)
        {
            cout << "out of bound" << endl;
            cout << "x: " << p0->X.transpose() << endl;
            cout << "xb: " << p0->Xb.transpose() << endl;
             cout << "xm: " << xm.transpose() << endl;
             p0->X = xm;
            int i = round(xm(0));
            int j = round(xm(1));
            int k = round(xm(2));
            for (size_t m=0; m<DomLBM->Flag[i][j][k].size(); ++m)
            {
                size_t ind = DomLBM->Flag[i][j][k][m];
                if (ind>=6){
                    double dis2 = (DomDEM->Lp[ind]->X-p0->X).squaredNorm();
                    double r = DomDEM->Lp[ind]->R;
                    double r2 = DomDEM->Lp[ind]->R*DomDEM->Lp[ind]->R;

                    Vector3d xc, n;
                    DomRWM->ContactPointWithSphere(p0, DomDEM->Lp[ind]->X, DomDEM->Lp[ind]->Xbr, DomDEM->Lp[ind]->R, xc, n);

                    cout << "ind " << ind << endl;
                    cout << DomDEM->Lp[ind]->R << endl;
                    cout << DomDEM->Lp[ind]->X.transpose() << endl;
                    cout << DomDEM->Lp[ind]->Xbr.transpose() << endl;

                    cout << "xc: " << xc.transpose() << endl;
                    cout << "n: " << n.transpose() << endl;
                    Vector3d v = -p0->Xb+xc + 2.*n.dot(p0->Xb-xc)*n;
                    cout << "v: " << v.transpose() << endl;
                    Vector3d x2 = (p0->X-xc).norm()*v.normalized()+xc;
                    cout << "x2: " << x2.transpose() << endl;
                }
            }
            abort();
        }
    }
    // for (size_t p=0; p<DomRWM->Lp.size(); ++p)
    // {
    //     RWM_PARTICLE* p0 = DomRWM->Lp[p];
    //     double dis = (p0->X-DomDEM->Lp[6]->X).norm()-DomDEM->Lp[6]->R;
    //     if (dis<0.)
    //     {
    //         cout << "still inside didnt pass test" << endl;
    //         cout << "p0->Xb: " << p0->Xb.transpose() << endl;
    //         cout << "p0->Vc: " << p0->Vc.transpose() << endl;
    //         cout << "p0->X: " << p0->X.transpose() << endl;
    //         double disb = (DomDEM->Lp[6]->Xbr-p0->Xb).norm()-DomDEM->Lp[6]->R;
    //          double disc = (DomDEM->Lp[6]->X-p0->Vc).norm()-DomDEM->Lp[6]->R;
    //         cout << "disb: " << disb << endl;
    //         cout << "dis: " << dis << endl;
    //         cout << "disc: " << disc << endl;

    //                 p0->X = p0->Vc;
    //                 Vector3d xc, n;
    //                 DomRWM->ContactPointWithSphere(p0, DomDEM->Lp[6]->X, DomDEM->Lp[6]->Xbr, DomDEM->Lp[6]->R, xc, n);
    //                 cout << "xc: " << xc.transpose() << endl;
    //                 cout << "n: " << n.transpose() << endl;
    //                 Vector3d v = -p0->Xb+xc + 2.*n.dot(p0->Xb-xc)*n;
    //                 cout << "v: " << v.transpose() << endl;
    //                 Vector3d x2 = (p0->X-xc).norm()*v.normalized()+xc;
    //                 cout << "x2: " << x2.transpose() << endl;
    //                 double dis2 = (x2-DomDEM->Lp[6]->X).norm()-DomDEM->Lp[6]->R;
    //                 cout << "dis2: " << dis2 << endl;
    //         abort();
    //     }
    // }
}

inline void DELBM::MoveVAM()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=6; i<DomDEM->Lp.size(); ++i)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[i];
        if (p0->Group==-1)
        {
            CrossPBC(p0);
            p0->VelocityVerlet(DomDEM->Dt);
            // if (updatebox)  p0->UpdateBox(D);
        }
        else
        {
            // #pragma omp critical
            for (size_t d=0; d<3; ++d)
            {
                #pragma omp atomic
                DomDEM->Lg[p0->Group]->Fh(d) += p0->Fh(d);
                #pragma omp atomic
                DomDEM->Lg[p0->Group]->Fc(d) += p0->Fc(d);
            }
        }            
    }
    for (size_t i=0; i<DomDEM->Lg.size(); ++i)
    {
        DEM_PARTICLE* g0 = DomDEM->Lg[i];
        CrossPBC(g0);
        g0->VelocityVerlet(DomDEM->Dt);
        // if (updatebox)  g0->UpdateBox(D);
        // cout << "g0->Fh= " << g0->Fh.transpose() << endl;
        // cout << "g0->Fc= " << g0->Fc.transpose() << endl;
        // cout << "g0->M*g0->G= " << g0->M*g0->G.transpose() << endl;
        // cout << "g0->V= " << g0->V.transpose() << endl;
        g0->ZeroForceTorque(true, true);
        Vector3d displace = g0->X - g0->Xb;
        for (size_t l=0; l<g0->Lp.size(); ++l)
        {
            size_t p = g0->Lp[l];
            DEM_PARTICLE* p0 = DomDEM->Lp[p];
            p0->X += displace;
            p0->V = g0->V;
            CrossPBC(p0);
            // if (updatebox)  p0->UpdateBox(D);
        }
    }
}

inline void DELBM::CrossPBC(DEM_PARTICLE* p0)
{
    bool cross = false;
    Vector3d x0 = p0->X;
    for (int d=0; d<DomLBM->D; ++d)
    {
        if (p0->X(d)>DomSize[d])
        {
            p0->X(d) = p0->X(d)-DomSize[d]-1;
            cross = true;
        }
        else if (p0->X(d)<0.)
        {
            p0->X(d) = p0->X(d)+DomSize[d]+1;
            cross = true;
        }
    }
    if (cross)
    {
        Vector3d add = p0->X-x0;
        for (size_t l=0; l<p0->Ln.size(); ++l)
        {
            p0->Ln[l][0] += (int) add(0);
            p0->Ln[l][1] += (int) add(1);
            p0->Ln[l][2] += (int) add(2);
        }
    }    
}

inline void DELBM::ApplyIBB(int method)
{
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        // cout << "p= " << p << endl;
        DEM_PARTICLE* p0 = DomDEM->Lp[p];

        for (size_t l=0; l<p0->Lb.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Lb[l][0];
            int j = DomDEM->Lp[p]->Lb[l][1];
            int k = DomDEM->Lp[p]->Lb[l][2];
            // cout << "i= " << i << " j= " << j << " k= " << k << endl;
            Vector3d ind (i,j,k);

            for (int q=0; q<DomLBM->Q; ++q)
            {
                // cout << "q= " << q << endl;
                Vector3d fh (0.,0.,0.);
                Vector3d th (0.,0.,0.);
                Vector3d n (0.,0.,0.);
                // cout << "aaaa" << endl;
                double delta = p0->Lq[l](q);
                // cout << " delta= " << delta << endl;
                if (delta>-1.)
                {
                    n = ind-delta*DomLBM->E[q]-p0->X;
                    Vector3d vp = p0->V + p0->W.cross(n);

                    if (method==0)          DomLBM->VIBB(i, j, k, q, delta, vp, fh);
                    else if (method==1)     DomLBM->IBBYu(i, j, k, q, delta, vp, fh);
                    
                    th = n.cross(fh);

                    p0->Fh += fh;
                    p0->Th += th;
                }
            } 
            // cout << "---------------" << endl;           
        }
    }
}

inline void DELBM::ApplyNEBB()
{
    // double Bexp = 1./1.5;
    double Bexp = 1.;
    // cout << "start NEBB" << endl; 
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        // cout << "p0->Lb.size()= " << p0->Lb.size() << endl;
        MatrixXd A(p0->Lb.size(), p0->Lb.size());
        A.setIdentity();
        MatrixXd B(p0->Lb.size(), 3);
        B.setZero();
        MatrixXd vwNorm(p0->Lb.size(), 3);
        vwNorm.setZero();
        MatrixXd X(p0->Lb.size(), 3);
        X.setZero();
        MatrixXd Arho(p0->Lb.size(), p0->Lb.size());
        Arho.setIdentity();
        MatrixXd Brho(p0->Lb.size(), 1);
        Brho.setZero();
        MatrixXd Xrho(p0->Lb.size(), 1);
        Xrho.setZero();

        vector<int> lbetween;
        lbetween.resize(0);

        for (size_t l=0; l<p0->Lb.size(); ++l)
        {
            int i = p0->Lb[l][0];
            int j = p0->Lb[l][1];
            int k = p0->Lb[l][2];

            Vector3d ind (i,j,k);

            double d = p0->Ld[l];
            Vector3d x = p0->Li[l];
            Vector3d n = (ind-p0->X).normalized();
            Vector3d vw = p0->V + p0->W.cross(p0->R*n);
            vwNorm.row(l) = vw.dot(n)*n;

            // Linear interpolation
            // Find node
            int xmin = int(x(0));
            int xmax = xmin+1;
            int ymin = int(x(1));
            int ymax = ymin+1;
            int zmin = int(x(2));
            int zmax = zmin+1;

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

            B.row(l) = (1.-pow(d/(d+Hi),Bexp))*(vw-vw.dot(n)*n);

            for (size_t ll=0; ll<pow(2,D); ++ll)
            {
                int in = (int) ver[ll](0);
                int jn = (int) ver[ll](1);
                int kn = (int) ver[ll](2);

                Vector3d s = x-ver[7-ll];
                double vol = abs(s(0)*s(1)*s(2));

                // Apply Periodic BC
                // if (DomLBM->Periodic[0])    in = (in+Nx+1)%(Nx+1);
                // if (DomLBM->Periodic[1])    jn = (jn+Ny+1)%(Ny+1);
                // if (DomLBM->Periodic[2])    kn = (kn+Nz+1)%(Nz+1);

                double gn = DomLBM->G[in][jn][kn][0];
                if ((gn-((int) gn))==0. && gn>-1.)
                {
                    lbetween.push_back(l);
                }
                if ((gn-((int) gn))==0.5 && (int) gn == p0->ID)
                {
                    vector< int > indn {in,jn,kn};
                    int ln = distance(p0->Lb.begin(), find (p0->Lb.begin(), p0->Lb.end(), indn));
                    A(l,ln) -= pow(d/(d+Hi),Bexp)*vol;
                    Arho(l,ln) -=vol;
                }
                else
                {
                    double rhon = 0.;
                    Vector3d vn (0., 0., 0.);
                    if ((gn-((int) gn))==0. && gn>-1.)
                    {
                        DEM_PARTICLE* p1 = DomDEM->Lp[(int) gn];
                        Vector3d xn (in, jn, kn);
                        vn = p1->V + p1->W.cross(xn-p1->X);
                        rhon = DomLBM->Rho0;
                    }
                    else
                    {
                        for (int q=0; q<DomLBM->Q; ++q)
                        {
                            rhon += DomLBM->F[in][jn][kn](q);
                            vn += DomLBM->F[in][jn][kn](q)*DomLBM->E[q];
                        }
                        vn /= rhon;
                    }

                    B.row(l) += pow(d/(d+Hi),Bexp)*vol*(vn-vn.dot(n)*n);
                    Brho(l,0) += vol*rhon;
                }
            }
        }
        // cout << "find martix" << endl; 
        // for (size_t lll=0; lll<lbetween.size(); ++lll)
        // {
        //     int m = lbetween(lll);
        //     A.row(m).setZero();
        //     A(m,m) = 1.;

        //     int i = p0->Lb[m][0];
        //     int j = p0->Lb[m][1];
        //     int k = p0->Lb[m][2];

        //     Vector3d ind (i,j,k);

        //     double d = p0->Ld[l];
        //     Vector3d x = p0->Li[l];
        //     Vector3d n = (ind-p0->X).normalized();
        //     Vector3d vw = p0->V + p0->W.cross(p0->R*n);
        //     vwNorm.row(l) = vw.dot(n)*n;
        //     B.row(m) = 
        // }
        // cout << "A= " << endl;
        // cout << A << endl;

        SparseLU<SparseMatrix<double>> dec0;
        dec0.compute(A.sparseView());
        X = dec0.solve(B); 

        // cout << "solve martix 11111111111" << endl; 

        SparseLU<SparseMatrix<double>> dec1;
        dec1.compute(Arho.sparseView());
        Xrho = dec1.solve(Brho);

        // cout << "solve martix 2222222222222" << endl; 

        for (size_t l=0; l<p0->Lb.size(); ++l)
        {
            int i = p0->Lb[l][0];
            int j = p0->Lb[l][1];
            int k = p0->Lb[l][2];

            Vector3d ind (i,j,k);

            double rhob = DomLBM->F[i][j][k](0);

            bool miss[DomLBM->Q];

            vector <Vector3d> lx;
            lx.resize(DomLBM->Q);

            vector <Vector3d> lvw;
            lvw.resize(DomLBM->Q);

            for (int q=0; q<DomLBM->Q; ++q)
            {
                miss[q] = false;
            }
            // bool between = false;
            // int idp1=-1;
            for (int q=1; q<DomLBM->Q; ++q)
            {
                // cout << "q= " << q << endl;
                int ip = (i - (int) DomLBM->E[q][0]);
                int jp = (j - (int) DomLBM->E[q][1]);
                int kp = (k - (int) DomLBM->E[q][2]);
                double gp = DomLBM->G[ip][jp][kp][0];
                // if (gp!=-1. && (gp-((int) gp))!=0.5)
                if ((int) gp==p0->ID)
                {
                    miss[q] = true;
                    // if (miss[q] && miss[DomLBM->Op[q]])
                    // {
                    //     // cout << "g= " << g << endl;
                    //     between = true;
                    //     if (idp1!=p0->ID)    idp1 = (int) gp;
                    // }
                }
            }

            rhob = Xrho(l,0);
            Vector3d vf = X.row(l)+vwNorm.row(l);
            VectorXd feq(DomLBM->Q);

            // if (between)
            // {
            //     // cout << "between" << endl;
            //     // cout << "idp1= " << idp1 << endl;
            //     DEM_PARTICLE* p1 = DomDEM->Lp[idp1];
            //     Vector3d n0 = (ind-p0->X).normalized();
            //     Vector3d vw0 = p0->V + p0->W.cross(p0->R*n0);
            //     double d0 = (ind-p0->X).norm()-p0->R;
            //     Vector3d n1 = (ind-p1->X).normalized();
            //     Vector3d vw1 = p1->V + p1->W.cross(p1->R*n1);
            //     double d1 = (ind-p1->X).norm()-p1->R;
            //     vf = d1/(d0+d1)*vw0 + d0/(d0+d1)*vw1;
            //     rhob = DomLBM->Rho0;
            //     DomLBM->CalFeqC(feq, rhob, vf);

            //     for (int q=1; q<DomLBM->Q; ++q)
            //     {
            //         if (miss[q])
            //         {
            //             // DomLBM->F[i][j][k](q) = DomLBM->Ft[i][j][k](DomLBM->Op[q]);
            //             DomLBM->F[i][j][k](q) = feq(q) + (DomLBM->Ft[i][j][k](DomLBM->Op[q]) - feq(DomLBM->Op[q]));
            //             // Vector3d fh = -(DomLBM->F[i][j][k](q) + DomLBM->Ft[i][j][k](DomLBM->Op[q]))*DomLBM->E[q];
            //             // // Vector3d th = (lx[q]-p0->X).cross(fh);
            //             // Vector3d th = (ind-p0->X).cross(fh);

            //             // p0->Fh += fh;
            //             // p0->Th += th;
            //         }
            //     }
            // }
            // else
            double g= DomLBM->G[i][j][k][0];
            if (g-((int) g)!=0.)
            {
                // cout << "into" << endl;
                DomLBM->CalFeqC(feq, rhob, vf);

                for (int q=1; q<DomLBM->Q; ++q)
                {
                    // int ip = (i - (int) DomLBM->E[q][0]);
                    // int jp = (j - (int) DomLBM->E[q][1]);
                    // int kp = (k - (int) DomLBM->E[q][2]);

                    // if (DomLBM->G[ip][jp][kp]==0.)
                    if (miss[q])
                    {
                        // cout << "miss" << endl;
                        DomLBM->F[i][j][k](q) = feq(q) + (DomLBM->F[i][j][k](DomLBM->Op[q]) - feq(DomLBM->Op[q]));

                        // F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);
                        // Vector3d fh = DomLBM->F[i][j][k](q)*(lvw[q]-DomLBM->E[q]) - DomLBM->Ft[i][j][k](DomLBM->Op[q])*(lvw[q]+DomLBM->E[q]);
                        Vector3d fh = -(DomLBM->F[i][j][k](q) + DomLBM->Ft[i][j][k](DomLBM->Op[q]))*DomLBM->E[q];
                        // Vector3d th = (lx[q]-p0->X).cross(fh);
                        // cout << "fh= " << fh.transpose() << endl;
                        // abort();
                        Vector3d th = (ind-p0->X).cross(fh);

                        p0->Fh += fh;
                        p0->Th += th;
                    }
                }
            }
            else if (p0->ID==0)
            {
            }
            // cout << "p0->Fh= " << p0->Fh.transpose() << endl;
        }
        // cout << "nebb" << endl; 
    }
}

inline void DELBM::ApplyLNEBB()
{
    for (size_t p=0; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];

        for (size_t l=0; l<p0->Lb.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Lb[l][0];
            int j = DomDEM->Lp[p]->Lb[l][1];
            int k = DomDEM->Lp[p]->Lb[l][2];

            Vector3d ind (i,j,k);

            double d = p0->Ld[l];
            Vector3d n = (ind-p0->X).normalized();
            Vector3d vw = p0->V + p0->W.cross(p0->R*n);
            Vector3d vwNorm = vw.dot(n)*n;
            Vector3d vwTang = vw-vwNorm;
            // Vector3d vfTang = (sqrt(3)-d)/sqrt(3)*vwTang + d/sqrt(3)*(DomLBM->V[i][j][k]-DomLBM->V[i][j][k].dot(n)*n);
            // Vector3d vfTang = (1.5-d)/1.5*vwTang + d/1.5*DomLBM->V[i][j][k];
            Vector3d vfTang = (1.5-d)*vwTang + (d-0.5)*(DomLBM->V[i][j][k]-DomLBM->V[i][j][k].dot(n)*n);
            Vector3d vf = vfTang+vwNorm;
            // Vector3d vf = vw;
            double rhof = DomLBM->Rho[i][j][k];

            VectorXd feq(DomLBM->Q);
            DomLBM->CalFeqC(feq, rhof, vf);

            for (int q=0; q<DomLBM->Q; ++q)
            {
                int ip = (i - (int) DomLBM->E[q][0]);
                int jp = (j - (int) DomLBM->E[q][1]);
                int kp = (k - (int) DomLBM->E[q][2]);
                double g = DomLBM->G[ip][jp][kp][0];
                if (g!=-1. && (g-((int) g))!=0.5)
                {
                    DomLBM->F[i][j][k](q) = feq(q) + (DomLBM->F[i][j][k](DomLBM->Op[q]) - feq(DomLBM->Op[q]));
                    Vector3d fh = -(DomLBM->F[i][j][k](q) + DomLBM->Ft[i][j][k](DomLBM->Op[q]))*DomLBM->E[q];
                    Vector3d th = (ind-p0->X).cross(fh);

                    p0->Fh += fh;
                    p0->Th += th;
                }
            }
        }
    } 
}

inline void DELBM::CalRhoVDELBM()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i = 0; i <= Nx; ++i)
    for (int j = 0; j <= Ny; ++j)
    for (int k = 0; k <= Nz; ++k)
    {
        DomLBM->Rho[i][j][k] = 0.;
        DomLBM->V[i][j][k] = Vector3d::Zero();

        for (int q = 0; q < DomLBM->Q; ++q)
        {
            DomLBM->Rho[i][j][k]    += DomLBM->F[i][j][k](q);
            DomLBM->V[i][j][k]      += DomLBM->F[i][j][k](q)*DomLBM->E[q];
        }

        DomLBM->V[i][j][k] /= DomLBM->Rho[i][j][k];
        double g = DomLBM->G[i][j][k][0];
        if (g>-1. && g-((int) g)==0.)
        {
            Vector3d ind (i,j,k);
            DEM_PARTICLE* p0 = DomDEM->Lp[(int) g];
            DomLBM->V[i][j][k] = p0->V + p0->W.cross(ind-p0->X);
        }
    }
}

inline void DELBM::UpdateXbrForRWM()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=6; p<DomDEM->Lp.size(); ++p)
    {
        DomDEM->Lp[p]->Xbr = DomDEM->Lp[p]->X;
    }
}

inline void DELBM::SolveOneStepVAM(int demNt, bool save, int ct)
{
    cout << "DomLBM->CollideSRT();" << endl;
    clock_t t_start = std::clock();
    DomLBM->CollideSRT();
    clock_t t_end = std::clock();
    cout << "CollideSRT time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    cout << "ApplyVAM;" << endl;
    t_start = std::clock();
    ApplyVAM();
    t_end = std::clock();
    cout << "ApplyVAM time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    cout << "DomLBM->Stream();" << endl;
    t_start = std::clock();
    DomLBM->Stream();
    t_end = std::clock();
    cout << "Stream time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // DomLBM->SetPeriodic(false, false, false);
    // DomLBM->SetWall();
    cout << "DomLBM->ApplyWall();" << endl;
    t_start = std::clock();
    DomLBM->ApplyWall();
    t_end = std::clock();
    cout << "ApplyWall time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    cout << "MoveVAM;" << endl;
    // DomDEM->Contact();
    t_start = std::clock();
    DomDEM->Dt = 1./demNt;
    bool saveFc = false;
    for (int demt=0; demt<demNt; ++demt)
    {
        DomDEM->UpdateXmirGlobal();
        if (demt==0)
        {
            UpdateXbrForRWM();
        }
        if (demt==demNt-1)
        {
            if (save)   saveFc = true;
        }
        DomDEM->Contact(saveFc, ct);
        // DomDEM->Move();
        MoveVAM();
        DomDEM->ZeroForceTorque(false, true);
    }
    DomDEM->Lc.clear();
    cout << "satert DomRWM->Move()" << endl;
    // if (UseRW)   DomRWM->Move();
    if (UseRW)   MoveRW();
    cout << "DomRWM->Move()" << endl;
    if (save)
    {
        // DomLBM->WriteFileH5(ct, 1);
        // DomDEM->WriteFileH5(ct);
        if (UseRW)  DomRWM->CalC();
        cout << "DomRWM->CalC()" << endl;
        WriteFileH5(ct, 1);
        cout << "WriteFileH5" << endl;
    }
    UpdateBoxGlobal();
    cout << "UpdateBoxGlobal" << endl;
    // DomDEM->WriteFileParticleInfo(ct);
    DomDEM->ZeroForceTorque(true, true);
    t_end = std::clock();
    cout << "MoveVAM time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;

    cout << "DomLBM->CalRhoV();" << endl;
    t_start = std::clock();
    DomLBM->CalRhoV();
    t_end = std::clock();
    cout << "CalRhoV time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
}

inline void DELBM::SolveVAM(int tt, int savet, int demNt)
{
    for (int t=0; t<tt; ++t)
    {
        bool save = false;
        if (t%savet==0)    save = true;
        SolveOneStepVAM(demNt, save, t);
    }
}

// inline void DELBM::SolveOneStepPSM(int demNt)
// {
//     DomLBM->CollideSRT();
//     ApplyPSM();
//     DomLBM->Stream();
//     DomLBM->SetPeriodic(false, false, false);
//     DomLBM->SetWall();
//     // DomDEM->Contact();
//     DomDEM->Dt = 1./demNt;
//     for (int demt=0; demt<demNt; ++demt)
//     {
//         // DomDEM->Contact();
//         // DomDEM->Move();
//         MoveVAM();
//         DomDEM->ZeroForceTorque(false, true);
//     }
//     DomDEM->ZeroForceTorque(true, true);
//     DomLBM->CalRhoV();
// }

inline void DELBM::SolveOneStepNEBB(int demNt)
{
    // cout << "DomLBM->CollideSRT();" << endl;
    // clock_t t_start = std::clock();
    DomLBM->CollideSRT();
    // clock_t t_end = std::clock();
    // cout << "CollideMRT time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomLBM->Stream();" << endl;
    // t_start = std::clock();
    DomLBM->Stream();
    // t_end = std::clock();
    // cout << "Stream time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomLBM->SetWall();" << endl;
    DomLBM->SetPeriodic(false, false, false);
    DomLBM->SetWall();
    // cout << "ApplyIBB(method);" << endl;
    // t_start = std::clock();
    ApplyNEBB();
    // t_end = std::clock();
    // cout << "ApplyNEBB time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomDEM->RecordX();" << endl;
    // DomDEM->RecordX();
    // cout << "DomDEM->Contact();" << endl;
    // DomDEM->Contact();
    // cout << "DomDEM->Move(1./demNt);" << endl;
    // t_start = std::clock();
    // DomDEM->Contact();

//     DEM_PARTICLE* p0 = DomDEM->Lp[0];
//     {
//         double delta = Nz-p0->x(2);
//     }

//     for (size_t i=0; i<DomDEM->Lc.size(); ++i)
//     {
//         DEM_PARTICLE* pi = DomDEM->Lp[DomDEM->Lc[i][0]];
//         DEM_PARTICLE* pj = DomDEM->Lp[DomDEM->Lc[i][1]];
//         DomDEM->Contact2P(pi,pj);
//         // Normal direction (pj pinnts to pi)
//         Vector3d n = pi->X-pj->X;
//         // Overlapping distance
//         double delta = n.norm()-pi->R-pj->R;
//         cout << "for lubrication delta= " << delta << endl;
//         double c = 0.5;
//         if (delta>0. && delta<c)
//         {
//                 double lamda = 1./(1./pi->R + 1./pj->R);
//                 n.normalize();
//                 Vector3d vr = (pi->V-pj->V).dot(n)*n;

//                 Vector3d fl = -6.*M_PI*DomLBM->Nu*DomLBM->Rho0*vr/(lamda*lamda)*(1./delta - 1./c);
//                 cout << "=====================" << endl;
//                 cout << "fl.transpose()= " << fl.transpose() << endl;
//                 cout << "=====================" << endl;
// /*                double lamda = 0.5*(1./pi->R + 1./pj->R);
//                 n.normalize();
//                 Vector3d vr = (pi->V-pj->V).dot(n)*n;
//                 Vector3d fl = -1.5*M_PI*DomLBM->Nu*DomLBM->Rho0*vr/(lamda*lamda*delta);*/
//                 pi->Fh += fl;
//                 pj->Fh -= fl;
//         }
//     }

    DomDEM->Dt = 1./demNt;
    for (int demt=0; demt<demNt; ++demt)
    {
        DomDEM->Move();
    }

    // cout << "pi->Fc.transpose()= " << DomDEM->Lp[0]->Fc.transpose() << endl;
    // cout << "pi->Fh.transpose()= " << DomDEM->Lp[0]->Fh.transpose() << endl;
    // t_end = std::clock();
    // cout << "Move time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "=================================" << endl;
    // cout << "fh= " << DomDEM->Lp[0]->Fh.transpose()/DomDEM->Lp[0]->M << endl;
    // cout << "fc= " << DomDEM->Lp[0]->Fc.transpose()/DomDEM->Lp[0]->M << endl;
    // cout << "=================================" << endl;

    // cout << "DomDEM->ZeroForceTorque();" << endl;
    DomDEM->ZeroForceTorque(true, true);
    // cout << "DomLBM->CalRhoV();" << endl;
    // t_start = std::clock();
    // CalRhoVDELBM();
    DomLBM->CalRhoV();
    // t_end = std::clock();
    // cout << "CalRhoV time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomLBM->V[25][25][20].norm()= " << DomLBM->V[25][25][20].norm() << endl;
    // cout << "UpdateG();" << endl;
    // t_start = std::clock();
    UpdateG();
    // t_end = std::clock();
    // cout << "UpdateG time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // Refill();
    // t_start = std::clock();
    // cout << "EqRefill();" << endl;
    EqRefill();
    // t_end = std::clock();
    // cout << "NormRefill time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "EqRefill();" << endl;
    // EqRefill();    
}

inline void DELBM::SolveOneStepIBB(int method, int demNt)
{
    // cout << "DomLBM->CollideMRT();" << endl;
    DomLBM->CollideSRT();
    // cout << "DomLBM->Stream();" << endl;
    DomLBM->Stream();
    // cout << "DomLBM->SetWall();" << endl;
    DomLBM->SetPeriodic(false, false, false);
    DomLBM->SetWall();
    // cout << "ApplyIBB(method);" << endl;
    ApplyIBB(method);
    // cout << "DomDEM->RecordX();" << endl;
    DomDEM->RecordX();
    // cout << "DomDEM->Contact();" << endl;
    // DomDEM->Contact();
    // cout << "DomDEM->Move(1./demNt);" << endl;
    DomDEM->Dt = 1./demNt;
    for (int demt=0; demt<demNt; ++demt)
    {
        DomDEM->Move();
    }
    // cout << "DomDEM->ZeroForceTorque();" << endl;
    DomDEM->ZeroForceTorque(true, true);
    // cout << "DomLBM->CalRhoV();" << endl;
    DomLBM->CalRhoV();
    // cout << "UpdateG();" << endl;
    UpdateGForIBB();
    // Refill();
    // NormRefill();
    // cout << "EqRefill();" << endl;
    EqRefill();
    // cout << "finish EqRefill();" << endl;
}

// inline void DELBM::WriteFileH5(int n, int scale)
// {
//     stringstream    out;                    //convert int to string for file name.
//     out << setw(6) << setfill('0') << n;            
//     string file_name_h5 = "DELBM"+out.str()+".h5";

//     H5File  file(file_name_h5, H5F_ACC_TRUNC);      //create a new hdf5 file.

//     hsize_t nx = (Nx+1)/scale;
//     hsize_t ny = (Ny+1)/scale;
//     hsize_t nz = (Nz+1)/scale;

//     int numLat = nx*ny*nz;

//     hsize_t dims_scalar[3] = {nx, ny, nz};          //create data space.
//     hsize_t dims_vector[4] = {nx, ny, nz, 3};       //create data space.

//     size_t npar = DomDEM->Lp.size()-6;
//     hsize_t dims_scalar_p[1] = {npar};           //create data space.
//     hsize_t dims_vector_p[1] = {3*npar};         //create data space.

//     int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
//     int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

//     int rank_scalar_p = sizeof(dims_scalar_p) / sizeof(hsize_t);
//     int rank_vector_p = sizeof(dims_vector_p) / sizeof(hsize_t);

//     DataSpace   *space_scalar = new DataSpace(rank_scalar, dims_scalar);
//     DataSpace   *space_vector = new DataSpace(rank_vector, dims_vector);

//     DataSpace   *space_scalar_p = new DataSpace(rank_scalar_p, dims_scalar_p);
//     DataSpace   *space_vector_p = new DataSpace(rank_vector_p, dims_vector_p);

//     double* rho_h5  = new double[numLat];
//     double* g_h5    = new double[numLat];
//     double* vel_h5  = new double[3*numLat];

//     double* r_p_h5    = new double[  npar];
//     double* rho_p_h5  = new double[  npar];
//     double* tag_p_h5  = new double[  npar];
//     double* pos_p_h5  = new double[3*npar];
//     double* vel_p_h5  = new double[3*npar];
//     double* agv_p_h5  = new double[3*npar];
//     double* fh_p_h5   = new double[3*npar];

//     int len = 0;
//     // #pragma omp parallel for schedule(static) num_threads(Nproc)
//     for (size_t k=0; k<nz; k++)
//     for (size_t j=0; j<ny; j++)
//     for (size_t i=0; i<nx; i++)
//     {       
//         rho_h5[len]     = 0.;
//         g_h5[len]       = 0.;
//         vel_h5[3*len  ] = 0.;
//         vel_h5[3*len+1] = 0.;
//         vel_h5[3*len+2] = 0.;
//         int cout = 0;
//         for (int kk=0; kk<scale; kk++)
//         for (int jj=0; jj<scale; jj++)
//         for (int ii=0; ii<scale; ii++)
//         {
//             if (i+ii<=(size_t) Nx && j+jj<=(size_t) Ny && k+kk<=(size_t) Nz)
//             {
//                 rho_h5[len]     += DomLBM->Rho[i+ii][j+jj][k+kk];
//                 g_h5[len]       += DomLBM->G[i+ii][j+jj][k+kk][0];
//                 vel_h5[3*len  ] += DomLBM->V[i][j][k](0);
//                 vel_h5[3*len+1] += DomLBM->V[i][j][k](1);
//                 vel_h5[3*len+2] += DomLBM->V[i][j][k](2);
//                 cout++;
//             }
//         }
//         rho_h5[len] /= (double) cout;
//         g_h5[len] /= (double) cout;
//         vel_h5[3*len  ] /= (double) cout;
//         vel_h5[3*len+1] /= (double) cout;
//         vel_h5[3*len+2] /= (double) cout;
//         len++;
//     }

//     for (size_t i=6; i<DomDEM->Lp.size(); ++i)
//     {
//         size_t j = i-6;
//         Vector3d agv = DomDEM->Lp[i]->Q._transformVector(DomDEM->Lp[i]->W);
//         r_p_h5  [  j  ]   = DomDEM->Lp[i]->R/(double) scale;
//         rho_p_h5[  j  ]   = DomDEM->Lp[i]->Rho;
//         tag_p_h5[  j  ]   = DomDEM->Lp[i]->Tag;
//         pos_p_h5[3*j  ]   = DomDEM->Lp[i]->X(0)/(double) scale;
//         pos_p_h5[3*j+1]   = DomDEM->Lp[i]->X(1)/(double) scale;
//         pos_p_h5[3*j+2]   = DomDEM->Lp[i]->X(2)/(double) scale;
//         vel_p_h5[3*j  ]   = DomDEM->Lp[i]->V(0);
//         vel_p_h5[3*j+1]   = DomDEM->Lp[i]->V(1);
//         vel_p_h5[3*j+2]   = DomDEM->Lp[i]->V(2);
//         fh_p_h5[3*j  ]    = DomDEM->Lp[i]->Fh(0);
//         fh_p_h5[3*j+1]    = DomDEM->Lp[i]->Fh(1);
//         fh_p_h5[3*j+2]    = DomDEM->Lp[i]->Fh(2);
//         agv_p_h5[3*j  ]   = agv(0);
//         agv_p_h5[3*j+1]   = agv(1);
//         agv_p_h5[3*j+2]   = agv(2);
//     }

//     DataSet *dataset_rho = new DataSet(file.createDataSet("Density", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet *dataset_g   = new DataSet(file.createDataSet("Gamma", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet *dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

//     DataSet *dataset_r_p      = new DataSet(file.createDataSet("Particle_Radius", PredType::NATIVE_DOUBLE, *space_scalar_p));
//     DataSet *dataset_rho_p    = new DataSet(file.createDataSet("Particle_Rho", PredType::NATIVE_DOUBLE, *space_scalar_p));
//     DataSet *dataset_tag_p    = new DataSet(file.createDataSet("Particle_Tag", PredType::NATIVE_DOUBLE, *space_scalar_p));
//     DataSet *dataset_pos_p    = new DataSet(file.createDataSet("Particle_Position", PredType::NATIVE_DOUBLE, *space_vector_p));
//     DataSet *dataset_vel_p    = new DataSet(file.createDataSet("Particle_Velocity", PredType::NATIVE_DOUBLE, *space_vector_p));
//     DataSet *dataset_agv_p    = new DataSet(file.createDataSet("Particle_AngularVelocity", PredType::NATIVE_DOUBLE, *space_vector_p));
//     DataSet *dataset_fh_p     = new DataSet(file.createDataSet("Particle_HydroForce", PredType::NATIVE_DOUBLE, *space_vector_p));

//     dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
//     dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
//     dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

//     dataset_r_p->write(r_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_rho_p->write(rho_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_tag_p->write(tag_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_pos_p->write(pos_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_vel_p->write(vel_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_agv_p->write(agv_p_h5, PredType::NATIVE_DOUBLE);
//     dataset_fh_p->write(fh_p_h5, PredType::NATIVE_DOUBLE);

//     delete space_scalar;
//     delete dataset_rho;
//     delete dataset_g;
//     delete dataset_vel;

//     delete rho_h5;
//     delete g_h5;
//     delete vel_h5;

//     // delete space_scalar;
//     // delete space_vector;

//     // delete dataset_r_p;
//     // delete dataset_rho_p;
//     // delete dataset_tag_p;
//     // delete dataset_pos_p;
//     // delete dataset_vel_p;
//     // delete dataset_agv_p;
//     // delete dataset_fh_p;

//     // delete r_p_h5;
//     // delete rho_p_h5;
//     // delete tag_p_h5;
//     // delete pos_p_h5;
//     // delete vel_p_h5;
//     // delete agv_p_h5;
//     // delete fh_p_h5;

//     file.close();

//     string file_name_xmf = "DELBM_"+out.str()+".xmf";

//     std::ofstream oss;
//     oss.open(file_name_xmf);
//     oss << "<?xml version=\"1.0\" ?>\n";
//     oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
//     oss << "<Xdmf Version=\"2.0\">\n";
//     oss << " <Domain>\n";
//     oss << "   <Grid Name=\"LBM\" GridType=\"Uniform\">\n";
//     oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << nz << " " << ny << " " << nx << "\"/>\n";
//     oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
//     oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
//     oss << "       </DataItem>\n";
//     oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << 1.0 << " " << 1.0  << " " << 1.0  << "\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Geometry>\n";
//     oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Density \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Velocity \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Gamma \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "   </Grid>\n";
//     oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
//     oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << DomDEM->Lp.size() << "\"/>\n";
//     oss << "     <Geometry GeometryType=\"XYZ\">\n";
//     oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << DomDEM->Lp.size() << " 3\" >\n";
//     oss << "        " << file_name_h5 <<":/Particle_Position \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Geometry>\n";
//     oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_Radius \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_Rho \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_Tag \n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_Velocity\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_AngularVelocity\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "     <Attribute Name=\"HydroForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
//     oss << "       <DataItem Dimensions=\"" << DomDEM->Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
//     oss << "        " << file_name_h5 <<":/Particle_HydroForce\n";
//     oss << "       </DataItem>\n";
//     oss << "     </Attribute>\n";
//     oss << "   </Grid>\n";
//     oss << " </Domain>\n";
//     oss << "</Xdmf>\n";
//     oss.close();
// }

inline void DELBM::WriteFileH5(int n, int scale)
{
    stringstream    out;                    //convert int to string for file name.
    out << setw(6) << setfill('0') << n;            
    string file_name_h5 = "DELBM"+out.str()+".h5";

    H5File  file(file_name_h5, H5F_ACC_TRUNC);      //create a new hdf5 file.

    hsize_t nx = (Nx+1)/scale;
    hsize_t ny = (Ny+1)/scale;
    hsize_t nz = (Nz+1)/scale;

    int numLat = nx*ny*nz;

    hsize_t dims_scalar[3] = {nx, ny, nz};          //create data space.
    hsize_t dims_vector[4] = {nx, ny, nz, 3};       //create data space.

    size_t npar = DomDEM->Lp.size()-6;

    hsize_t dims_scalar_p[1] = {npar};            //create data space.
    hsize_t dims_vector_p[1] = {3*npar};          //create data space.

    int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
    int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

    int rank_scalar_p = sizeof(dims_scalar_p) / sizeof(hsize_t);
    int rank_vector_p = sizeof(dims_vector_p) / sizeof(hsize_t);

    DataSpace   *space_scalar = new DataSpace(rank_scalar, dims_scalar);
    DataSpace   *space_vector = new DataSpace(rank_vector, dims_vector);

    DataSpace   *space_scalar_p = new DataSpace(rank_scalar_p, dims_scalar_p);
    DataSpace   *space_vector_p = new DataSpace(rank_vector_p, dims_vector_p);

    double* rho_h5  = new double[numLat];
    double* g_h5    = new double[numLat];
    double* vel_h5  = new double[3*numLat];

    double* r_h5_p    = new double[  npar];
    double* rho_h5_p  = new double[  npar];
    double* tag_h5_p  = new double[  npar];
    double* pos_h5_p  = new double[3*npar];
    double* vel_h5_p  = new double[3*npar];
    double* agv_h5_p  = new double[3*npar];
    double* fh_h5_p   = new double[3*npar];

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
                rho_h5[len] += DomLBM->Rho[i+ii][j+jj][k+kk];
                g_h5[len] += DomLBM->G[i+ii][j+jj][k+kk][0];
                vel_h5[3*len  ] += DomLBM->V[i][j][k](0);
                vel_h5[3*len+1] += DomLBM->V[i][j][k](1);
                vel_h5[3*len+2] += DomLBM->V[i][j][k](2);
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

    for (size_t i=0; i<npar; ++i)
    {
        size_t ind = i+6;
        Vector3d agv = DomDEM->Lp[ind]->Q._transformVector(DomDEM->Lp[ind]->W);
        r_h5_p  [  i  ]   = DomDEM->Lp[ind]->R;
        rho_h5_p[  i  ]   = DomDEM->Lp[ind]->Rho;
        tag_h5_p[  i  ]   = DomDEM->Lp[ind]->Tag;
        pos_h5_p[3*i  ]   = DomDEM->Lp[ind]->X(0);
        pos_h5_p[3*i+1]   = DomDEM->Lp[ind]->X(1);
        pos_h5_p[3*i+2]   = DomDEM->Lp[ind]->X(2);
        vel_h5_p[3*i  ]   = DomDEM->Lp[ind]->V(0);
        vel_h5_p[3*i+1]   = DomDEM->Lp[ind]->V(1);
        vel_h5_p[3*i+2]   = DomDEM->Lp[ind]->V(2);
        fh_h5_p[3*i  ]    = DomDEM->Lp[ind]->Fh(0);
        fh_h5_p[3*i+1]    = DomDEM->Lp[ind]->Fh(1);
        fh_h5_p[3*i+2]    = DomDEM->Lp[ind]->Fh(2);
        agv_h5_p[3*i  ]   = agv(0);
        agv_h5_p[3*i+1]   = agv(1);
        agv_h5_p[3*i+2]   = agv(2);
    }

    DataSet *dataset_rho = new DataSet(file.createDataSet("Density", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet *dataset_g   = new DataSet(file.createDataSet("Gamma", PredType::NATIVE_DOUBLE, *space_scalar));
    DataSet *dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_DOUBLE, *space_vector));

    DataSet *dataset_r_p      = new DataSet(file.createDataSet("Radius_P", PredType::NATIVE_DOUBLE, *space_scalar_p));
    DataSet *dataset_rho_p    = new DataSet(file.createDataSet("Rho_P", PredType::NATIVE_DOUBLE, *space_scalar_p));
    DataSet *dataset_tag_p    = new DataSet(file.createDataSet("Tag_P", PredType::NATIVE_DOUBLE, *space_scalar_p));
    DataSet *dataset_pos_p    = new DataSet(file.createDataSet("Position_P", PredType::NATIVE_DOUBLE, *space_vector_p));
    DataSet *dataset_vel_p    = new DataSet(file.createDataSet("Velocity_P", PredType::NATIVE_DOUBLE, *space_vector_p));
    DataSet *dataset_agv_p    = new DataSet(file.createDataSet("AngularVelocity_P", PredType::NATIVE_DOUBLE, *space_vector_p));
    DataSet *dataset_fh_p     = new DataSet(file.createDataSet("HydroForce_P", PredType::NATIVE_DOUBLE, *space_vector_p));

    dataset_rho->write(rho_h5, PredType::NATIVE_DOUBLE);
    dataset_g->write(g_h5, PredType::NATIVE_DOUBLE);
    dataset_vel->write(vel_h5, PredType::NATIVE_DOUBLE);

    dataset_r_p->write(r_h5_p, PredType::NATIVE_DOUBLE);
    dataset_rho_p->write(rho_h5_p, PredType::NATIVE_DOUBLE);
    dataset_tag_p->write(tag_h5_p, PredType::NATIVE_DOUBLE);
    dataset_pos_p->write(pos_h5_p, PredType::NATIVE_DOUBLE);
    dataset_vel_p->write(vel_h5_p, PredType::NATIVE_DOUBLE);
    dataset_agv_p->write(agv_h5_p, PredType::NATIVE_DOUBLE);
    dataset_fh_p->write(fh_h5_p, PredType::NATIVE_DOUBLE);

    if (UseRW)
    {
        double* c_h5    = new double[3*numLat];

        int len = 0;
        // #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t k=0; k<nz; k++)
        for (size_t j=0; j<ny; j++)
        for (size_t i=0; i<nx; i++)
        {       
            c_h5[len] = 0.;
            int cout = 0;
            for (int kk=0; kk<scale; kk++)
            for (int jj=0; jj<scale; jj++)
            for (int ii=0; ii<scale; ii++)
            {
                if (i+ii<=(size_t) Nx && j+jj<=(size_t) Ny && k+kk<=(size_t) Nz)
                {
                    c_h5[len] += DomRWM->C[i+ii][j+jj][k+kk];
                    cout++;
                }
            }
            c_h5[len] /= (double) cout;
            len++;
        }
        DataSet *dataset_c = new DataSet(file.createDataSet("Concentration", PredType::NATIVE_DOUBLE, *space_scalar));
        dataset_c->write(c_h5, PredType::NATIVE_DOUBLE);
        delete dataset_c;
        delete c_h5;
    }

    delete space_scalar;
    delete dataset_rho;
    delete dataset_g;
    delete dataset_vel;
    delete rho_h5;
    delete g_h5;
    delete vel_h5;

    delete space_scalar_p;
    delete space_vector_p;
    delete dataset_r_p;
    delete dataset_rho_p;
    delete dataset_tag_p;
    delete dataset_pos_p;
    delete dataset_vel_p;
    delete dataset_agv_p;
    delete dataset_fh_p;

    file.close();

    string file_name_xmf = "DELBM_"+out.str()+".xmf";

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
    if (UseRW)
    {
        oss << "     <Attribute Name=\"Concentration\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << file_name_h5 <<":/Concentration \n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
    }
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
    oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << npar << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << npar << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position_P \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius_P\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Radius_P \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Rho_P\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Rho_P \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag_P\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Tag_P \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity_P\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity_P\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngularVelocity_P\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/AngularVelocity_P\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"HydroForce_P\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/HydroForce_P\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}