#include <LBM.h>

int main(int argc, char const *argv[])
{
	int nx = 50;
	int ny = 50;
	int nz = 0;

	double nu = 0.1;

	Vector3d ac (1.0e-6, 0., 0.);

/*===============================================================================*/

	LBM* lbm = new LBM(D2Q9, SRT, false, nx ,ny ,nz , nu);

	Vector3d v0 (0.,0.,0.);
	lbm->Init(1., v0);

	lbm->SetA(ac);

	for (size_t i=0; i<=nx; ++i)
	{
		Vector3i x_top (i,ny,0);
		Vector3i x_bot (i,0 ,0);
		lbm->Lwall.push_back(x_top);
		lbm->Lwall.push_back(x_bot);
	}

	for (int t=0; t<10000; ++t)
	{
		if (t%100==0)
		{
			cout << "Time= " << t << endl;
			lbm->WriteFileH5(t,1);
		}
		lbm->CollideSRT();
		lbm->ApplyWall();
		lbm->Stream();
		lbm->CalRhoV();	
	}

	return 0;
}