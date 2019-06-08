#include "../HEADER.h"
#include <DEM.h>
#include <LBM.h>

class DELBM
{
public:
	DELBM();
	~DELBM();
	DELBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu);
	void Init(double rho0, Vector3d initV);
	void InitG(DEM_PARTICLE* p0);
	void UpdateG();
	void UpdateLn(DEM_PARTICLE* p0);
    void Refill();
    void EqRefill();
    void NormRefill();
    void ApplyIBB(int method);                                                          // "0" for VIBB, "1" for Yu's double IBB
    void ApplyNEBB();
    void ApplyLNEBB();
    double CalCrossPointQ(Vector3d& a, Vector3d& v, Vector3d& c, double r);
	void SloveOneStep(int method, int demNt);
    void SloveOneStepNEBB(int demNt);
    void CalRhoVDELBM();

	LBM* 							DomLBM;												// Domain of LBM
	DEM*							DomDEM;												// Domain of DEM

    int 							Nx;													// Domain size
    int 							Ny;
    int 							Nz;

	int 							Nproc;												// Number of processors which used
    int 							D;

    double                          Hi;                                                 // Interpolation distance
    int                             BCmethod;                                           // Boundary condition methods, 0 for NEBB and 1 for IMB
    int                             Gmethod;                                            // Methods for calculate Gamma for IMB, 0 for distance based method and 1 for subgrid method

    vector< vector<int> >			Lr;                                                 // List of refilling nodes
};

inline DELBM::DELBM(DnQm dnqm, CollisionModel cmodel, bool incompressible, int nx, int ny, int nz, double nu)
{
	DomLBM = new LBM(dnqm, cmodel, incompressible, nx, ny, nz, nu);
	DomDEM = new DEM(nx, ny, nz);

	Nx = nx;
	Ny = ny;
	Nz = nz;

    D = DomLBM->D;

    Hi = 0.25;

    Nproc = 12;
}

inline void DELBM::Init(double rho0, Vector3d initV)
{
	DomLBM->Init(rho0, initV);

	for (size_t p=0; p<DomDEM->Lp.size(); ++p)
    {
        InitG(DomDEM->Lp[p]);
    }

    DomLBM->Nproc = Nproc;

    Lr.resize(0);
    DomDEM->Lc.resize(0);
}

inline void DELBM::InitG(DEM_PARTICLE* p0)
{
    int xMin, xMax, yMin, yMax, zMin, zMax;
    // Only init Gamma inside of the surrounding box for particles.
    xMin = p0->MinX;
    xMax = p0->MaxX;

    yMin = p0->MinY;
    yMax = p0->MaxY;

    zMin = p0->MinZ;
    zMax = p0->MaxZ;
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
                    int p = p0->ID;
                    int q = (int) DomLBM->G[i][j][k][0];
                    DomDEM->Lc.push_back({min(p,q), max(p,q)});
                }
            }
        }
    }
    // For IMB
    else if (BCmethod==1)
    {
        for (int i=xMin; i<=xMax; ++i)
        for (int j=yMin; j<=yMax; ++j)
        for (int k=zMin; k<=zMax; ++k)
        {
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
                    int p0 = (int) G[i][j][k][0];
                    int p1 = p0->ID;
                    Lc.push_back({min(p0,p1), max(p0,p1)});
                }
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
}

inline void DELBM::UpdateG()
{
    DomDEM->Lc.resize(0);
    // Reset flag for neighbor nodes.
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t p=0; p<DomDEM->Lp.size(); ++p)
    {
        for (size_t l=0; l<DomDEM->Lp[p]->Ln.size(); ++l)
        {
            int i = DomDEM->Lp[p]->Ln[l][0];
            int j = DomDEM->Lp[p]->Ln[l][1];
            int k = DomDEM->Lp[p]->Ln[l][2];
            // Reset neighour nodes to fluid nodes (g=-1)
            if (i>=0 && i<=Nx && j>=0 && j<=Ny && k>=0 && k<=Nz)
            {
                if (DomLBM->G[i][j][k][0]==-2.)
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
                else
                {
                    if ((g-((int) g)) == 0.5 && (id-((int) id))==0)
                    {
                        DomLBM->G[i][j][k][0] = id;
                    }

                    int p = p0->ID;
                    int q = (int) DomLBM->G[i][j][k][0];
                    if (q==-2)
                    {
                        cout << "Contact" << endl;
                        Vector3d n = ind-DomDEM->Lp[p]->X;
                        // Overlapping distance
                        double delta = DomDEM->Lp[p]->R+0.5-n.norm();
                        if (delta>0)
                        {
                            n.normalize();
                            Vector3d fc = n*pow(delta, 1.5);
                            DomDEM->Lp[p]->Fc += DomDEM->Lp[p]->Kn*fc;
                        }
                    }
                    else
                    {
                        #pragma omp critical
                        {
                            if (p!=q)   DomDEM->Lc.push_back({min(p,q), max(p,q)});
                        }
                    }
                }
            }
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

inline void DELBM::UpdateG2()
{
    // Reset flag for neighbor nodes.
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int p=0; p<DomDEM->Lp.size(); p++)
    {
         for (int l=0; l<DomDEM->Lp[p]->Ln.size(); l++)
        {
            int i = DomDEM->Lp[p]->Ln[l][0];
            int j = DomDEM->Lp[p]->Ln[l][1];
            int k = DomDEM->Lp[p]->Ln[l][2];

            if (DomLBM->G[i][j][k][0] > -2.)
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
    for (int p=0; p<DomDEM->Lp.size(); p++)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];
        // Reset boundary nodes list.
        p0->Lb.clear();
        for (int l=0; l<p0->Ln.size(); l++)
        {
            int i = p0->Ln[l][0];
            int j = p0->Ln[l][1];
            int k = p0->Ln[l][2];

            double tempG = 0;
            // if (DomDEM->Lp[p]->Type==0)  tempG = CalSphereG(i,j,k,DomDEM->Lp[p]);
            if (Gmethod==0)            tempG = CalGdistance(i,j,k,p0);
            else if (Gmethod==1)       tempG = CalGvolumeSphere(i,j,k,p0, 20);
            else
            {
                cout << "wrong Gmethod!" << endl;
                abort();
            }
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

            if (i>=0 && i<=Nx && j>=0 && j<=Ny && k>=0 && k<=Nz)
            {
                if (DomLBM->G[i][j][k][0]!=-2.)
                {
                    #pragma omp critical
                    {
                        p0->Ln.push_back({i,j,k});
                    }
                }
            }
        }
    }
    sort( p0->Ln.begin(), p0->Ln.end() );
    p0->Ln.erase( unique( p0->Ln.begin(), p0->Ln.end() ), p0->Ln.end() );
}

inline double DELBM::CalGdistance(int i, int j, int k, DEM_PARTICLE* p0)
{
    Vector3d xc (i,j,k);
    double dis = abs((xc-p0->X).norm()-p0->R);
    double g = 0.;

    if (dis<1.)         g = 1.-dis;
    // else if (dis<1.) g = 1.-dis;
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

        Vector3d ind (i,j,k);

        int p = (int) DomLBM->G[i][j][k][0];

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

            if (DomLBM->G[ip][jp][kp][0]==0.)
            {
                DomLBM->F[i][j][k](q) = feq(q);
            }
        }
    }
}

inline void DELBM::ApplyIBB(int method)
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

            for (int q=0; q<DomLBM->Q; ++q)
            {
                Vector3d fh (0.,0.,0.);
                Vector3d th (0.,0.,0.);
                Vector3d n (0.,0.,0.);
                double delta = p0->Lq[l](q);

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
        }
    }
}

inline void DELBM::ApplyNEBB()
{
    double Bexp = 1./1.5;
    // double Bexp = 1.;
    // cout << "start NEBB" << endl; 
    for (size_t p=0; p<DomDEM->Lp.size(); ++p)
    {
        DEM_PARTICLE* p0 = DomDEM->Lp[p];

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

        SparseLU<SparseMatrix<double>> dec0;
        dec0.compute(A.sparseView());
        X = dec0.solve(B); 

        SparseLU<SparseMatrix<double>> dec1;
        dec1.compute(Arho.sparseView());
        Xrho = dec1.solve(Brho);

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
            bool between = false;
            int idp1=-1;
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
                DomLBM->CalFeqC(feq, rhob, vf);

                for (int q=1; q<DomLBM->Q; ++q)
                {
                    // int ip = (i - (int) DomLBM->E[q][0]);
                    // int jp = (j - (int) DomLBM->E[q][1]);
                    // int kp = (k - (int) DomLBM->E[q][2]);

                    // if (DomLBM->G[ip][jp][kp]==0.)
                    if (miss[q])
                    {
                        DomLBM->F[i][j][k](q) = feq(q) + (DomLBM->F[i][j][k](DomLBM->Op[q]) - feq(DomLBM->Op[q]));

                        // F[i][j][k](q)*(vw-E[q]) - Ft[i][j][k](Op[q])*(E[q]+vw);
                        // Vector3d fh = DomLBM->F[i][j][k](q)*(lvw[q]-DomLBM->E[q]) - DomLBM->Ft[i][j][k](DomLBM->Op[q])*(lvw[q]+DomLBM->E[q]);
                        Vector3d fh = -(DomLBM->F[i][j][k](q) + DomLBM->Ft[i][j][k](DomLBM->Op[q]))*DomLBM->E[q];
                        // Vector3d th = (lx[q]-p0->X).cross(fh);
                        Vector3d th = (ind-p0->X).cross(fh);

                        p0->Fh += fh;
                        p0->Th += th;
                    }
                }
            }
            else if (p0->ID==0)
            {
/*                DEM_PARTICLE* p1 = DomDEM->Lp[(int) g];
                double lamda = 0.5*(1./p0->R + 1./p1->R);
                double delta = (p0->X-p1->X).norm()-p0->R-p1->R;
                Vector3d n = (p1->X-p0->X).normalized();
                Vector3d vr = (p1->V-p0->V).dot(n)*n;
                double c = 1.;
                Vector3d fl = 1.5*M_PI*DomLBM->Nu*DomLBM->Rho0*vr/(lamda*lamda)*(1./delta-1./c);
                cout << "=====================" << endl;
                cout << "fl.transpose()= " << fl.transpose() << endl;
                cout << "=====================" << endl;
                // double fl = 1.5*0.6*DomLBM->Nu*DomLBM->Rho0*vr.norm()/lamda/(1./delta/delta - 1./c/c)/(c*c*c);

                // p0->Fh += fl;
                for (int q=1; q<DomLBM->Q; ++q)
                {
                    if (miss[q])
                    {
                        p0->Fh += fl;
                        // p0->Fh += fl*DomLBM->E[q];
                            // p0->Th += th;
                    }
                }*/
                
                // Vector3d n = (ind-p0->X).normalized();
                // Vector3d vw = p0->V + p0->W.cross(p0->R*n);
                // DEM_PARTICLE* p1 = DomDEM->Lp[(int) g];
                // Vector3d xt = p0->X + (p0->R+0.0001)*n;

                // if ((xt-p1->X).norm()>p1->R)
                // {
                //     vf = vw;
                //     rhob = DomLBM->Rho0;
                //     DomLBM->CalFeqC(feq, rhob, vf);
                //     for (int q=1; q<DomLBM->Q; ++q)
                //     {
                //         if (miss[q])
                //         {
                //             Vector3d fh = -DomLBM->Tau*(2.*feq(DomLBM->Op[q]) + 6.*DomLBM->W[q]*rhob*DomLBM->E[q].dot(vw))*DomLBM->E[q];
                //             Vector3d th = (ind-p0->X).cross(fh);

                //             p0->Fh += fh;
                //             p0->Th += th;
                //         }
                //     }
                // }
            }
        }
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

inline void DELBM::SloveOneStepNEBB(int demNt)
{
    // cout << "DomLBM->CollideSRT();" << endl;
    // clock_t t_start = std::clock();
    DomLBM->CollideMRT();
    // clock_t t_end = std::clock();
    // cout << "CollideMRT time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomLBM->Stream();" << endl;
    // t_start = std::clock();
    DomLBM->Stream();
    // t_end = std::clock();
    // cout << "Stream time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "DomLBM->SetWall();" << endl;
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
    DomDEM->ZeroForceTorque();
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
    EqRefill();
    // t_end = std::clock();
    // cout << "NormRefill time= " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << endl;
    // cout << "EqRefill();" << endl;
    // EqRefill();    
}

inline void DELBM::SloveOneStep(int method, int demNt)
{
    // cout << "DomLBM->CollideSRT();" << endl;
/*    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i=0; i<=Nx; i++)
    for (int j=0; j<=Ny; j++)
    for (int k=0; k<=Nz; k++)
    {
        // cout << "i= " << i << " j= " << j << " k= " << k << endl;
        double g = DomLBM->G[i][j][k];
        // cout << "g= " << g << endl;

        if (g!=-1. && (g-((int) g)) != 0.5)
        {
            int p = (int) g;
            // cout << "p= " << p << endl;
            DEM_PARTICLE* p0 = DomDEM->Lp[p];

            Vector3d ind (i,j,k);

            // Vector3d v = p0->V + p0->W.cross(ind-p0->X);

            // cout << "v= " << v.transpose() << endl;
            // cout << "DomLBM->V[i][j][k]= " << DomLBM->V[i][j][k].transpose() << endl;

            DomLBM->V[i][j][k] = p0->V + p0->W.cross(ind-p0->X);

            // cout << "after    DomLBM->V[i][j][k]= " << DomLBM->V[i][j][k].transpose() << endl;
        }

        // cout << "DomLBM->CollideSRTLocal(i, j, k);" << endl;
        DomLBM->CollideSRTLocal(i, j, k);
        // cout << "DomLBM->ExForce[i][j][k]" << endl;
        DomLBM->ExForce[i][j][k] = Vector3d::Zero();
    }*/

    cout << "DomLBM->CollideMRT();" << endl;
    DomLBM->CollideMRT();
    cout << "DomLBM->Stream();" << endl;
    DomLBM->Stream();
    cout << "DomLBM->SetWall();" << endl;
    DomLBM->SetWall();
    cout << "ApplyIBB(method);" << endl;
    ApplyIBB(method);
    cout << "DomDEM->RecordX();" << endl;
    DomDEM->RecordX();
    cout << "DomDEM->Contact();" << endl;
    DomDEM->Contact();
    cout << "DomDEM->Move(1./demNt);" << endl;

    DomDEM->Dt = 1./demNt;
    for (int demt=0; demt<demNt; ++demt)
    {
        DomDEM->Move();
    }

    // cout << "=================================" << endl;
    // cout << "fh= " << DomDEM->Lp[0]->Fh.transpose()/DomDEM->Lp[0]->M << endl;
    // cout << "fc= " << DomDEM->Lp[0]->Fc.transpose()/DomDEM->Lp[0]->M << endl;
    // cout << "=================================" << endl;

    cout << "DomDEM->ZeroForceTorque();" << endl;
    DomDEM->ZeroForceTorque();
    cout << "DomLBM->CalRhoV();" << endl;
    DomLBM->CalRhoV();
    cout << "UpdateG();" << endl;
    UpdateG();
    // Refill();
    // NormRefill();
    cout << "EqRefill();" << endl;
    EqRefill();
}