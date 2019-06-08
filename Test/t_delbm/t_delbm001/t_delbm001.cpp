#include <DELBM.h>

int main(int argc, char const *argv[])
{
	double dx = 1.0e-3; 				/*m*/
	double dt = 2.0e-4; 				/*s*/

	double nx_physical = 0.1;			/*m*/
	double ny_physical = 0.1;			/*m*/
	double nz_physical = 0.16;			/*m*/

	double rhof_physical = 960.;		/*kg/m^3*/
	double rhos_physical = 1120.;		/*kg/m^3*/

	double r_physical = 0.0075;			/*m*/

	double miu_physical = 0.058;		/*N.s/m^2*/

	double nu_physical = miu_physical/rhof_physical;		/*m^2/s*/

	Vector3d pos_physical (0.5*nx_physical, 0.5*ny_physical, 0.123);

	Vector3d g_physical (0., 0., -9.8);

	// To make sure rhof=1.
	double dm = rhof_physical*dx*dx*dx; /*kg*/

/*===============================================================================*/

	int nx = (int) (nx_physical/dx);
	int ny = (int) (ny_physical/dx);
	int nz = (int) (nz_physical/dx);

	double rhof = rhof_physical*dx*dx*dx/dm;
	double rhos = rhos_physical*dx*dx*dx/dm;

	double r = r_physical/dx;
	double nu = nu_physical*dt/(dx*dx);
	Vector3d pos0 = pos_physical/dx;
	Vector3d g = (1.-rhof/rhos)*g_physical*dt*dt/dx;

	cout << " rhof= " << rhof << " rhos= " << rhos << endl;
	cout << "nx= " << nx << " ny= " << ny << " nz= " << nz << endl;
	cout << "nu= " << nu << endl;
	cout << "g= " << g.transpose() << endl;

/*===============================================================================*/

	DELBM* delbm = new DELBM(D3Q15, MRT, true, nx ,ny ,nz , nu);

	delbm->DomDEM->AddSphere(-1, r, pos0, rhos);

	Vector3d v0 (0.,0.,0.);
	delbm->Init(1., v0);

	delbm->DomLBM->SetPeriodic(false, false, false);

	delbm->DomDEM->SetG(g);

	for (int t=0; t<1000000; ++t)
	{
		// delbm->UpdateG();

		if (t%100==0)
		{
			cout << "Time= " << t << endl;
			delbm->DomLBM->WriteFileH5(t);
			delbm->DomDEM->WriteFileParticleInfo(t);
			cout << "x= " << delbm->DomDEM->Lp[0]->X.transpose() << endl;
			cout << "vz= " << delbm->DomDEM->Lp[0]->V(2)*dx/dt << endl;
		}	
		delbm->SloveOneStepNEBB(1);
	}

	return 0;
}