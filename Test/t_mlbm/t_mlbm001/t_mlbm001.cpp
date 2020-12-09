#include <MLBM.h>

int main(int argc, char const *argv[])
{
	int nx = 120;
	int ny = 90;
	int nz = 0;

	double nu = 0.1;
	double dc = 0.01;

	Vector3d ac (1.0e-6, 0., 0.);

	VectorXd pm(10);

	pm(0) = 0.;
	pm(1) = 1.;
	pm(2) = 2.;

	cout << pm.transpose() << endl;
	pm.resize(3);
	cout << pm.size() << endl;
	cout << pm.transpose() << endl;
	abort();
/*===============================================================================*/

	MLBM* mlbm = new MLBM(D2Q9, SRT, false, nx ,ny ,nz , nu, dc);

	Vector3d v0 (0.,0.,0.);
	mlbm->Init(1., 0., v0);

	cout << "finish init" << endl;

	// mlbm->DomF->SetA(ac);

	Vector3d x0 (15, 0.75*ny, 0.);
	double r = 5.;
	double c0= 1.;

	Vector3d vin (0.01, 0., 0.);

	// add concentration
	for (int i=0; i<=nx; ++i)
	for (int j=0; j<=ny; ++j)
	{
		Vector3d xl (i, j, 0);
		if ((xl-x0).norm()<r)
		{
			// mlbm->DomS->AddCLocal(i,j,0,c0);
			mlbm->DomS->AddSourceLocal(i,j,0,0.01*c0);
		}
	}

	int length = 50;

	// set velocity inlet
	for (int i=0; i<=0; ++i)
	for (int j=length+1; j<=ny-1; ++j)
	{
		VelBC nodeVel;
		nodeVel.Pos << i,j,0;
		nodeVel.Nei << i+1,j,0;
		nodeVel.Vel = vin;

		mlbm->DomF->Lvel.push_back(nodeVel);
	}

	// set non gradient outlet
	for (int i=nx; i<=nx; ++i)
	for (int j=length+1; j<=ny-1; ++j)
	{
		NoGradBC nodeNog;
		nodeNog.Pos << i,j,0;
		nodeNog.Nei << i-1,j,0;
		mlbm->DomF->Lnog.push_back(nodeNog);
		mlbm->DomS->Lnog.push_back(nodeNog);
	}

	// set wall
	for (int i=0; i<=1; ++i)
	for (int j=1; j<=length; ++j)
	{
		Vector3i x_b (i,j,0);
		mlbm->DomF->Lwall.push_back(x_b);
		mlbm->DomS->Lwall.push_back(x_b);
	}

	for (int i=29; i<=30; ++i)
	for (int j=ny-length; j<ny; ++j)
	{
		Vector3i x_b (i,j,0);
		mlbm->DomF->Lwall.push_back(x_b);
		mlbm->DomS->Lwall.push_back(x_b);
	}

	for (int i=59; i<=60; ++i)
	for (int j=1; j<=length; ++j)
	{
		Vector3i x_b (i,j,0);
		mlbm->DomF->Lwall.push_back(x_b);
		mlbm->DomS->Lwall.push_back(x_b);
	}

	for (int i=89; i<=90; ++i)
	for (int j=ny-length; j<ny; ++j)
	{
		Vector3i x_b (i,j,0);
		mlbm->DomF->Lwall.push_back(x_b);
		mlbm->DomS->Lwall.push_back(x_b);
	}

	for (int i=nx-1; i<=nx; ++i)
	for (int j=1; j<=length; ++j)
	{
		Vector3i x_b (i,j,0);
		mlbm->DomF->Lwall.push_back(x_b);
		mlbm->DomS->Lwall.push_back(x_b);
	}

	for (int i=0; i<=nx; ++i)
	{
		Vector3i x_top (i,ny,0);
		Vector3i x_bot (i,0 ,0);
		mlbm->DomF->Lwall.push_back(x_top);
		mlbm->DomF->Lwall.push_back(x_bot);
		mlbm->DomS->Lwall.push_back(x_top);
		mlbm->DomS->Lwall.push_back(x_bot);
	}
	cout << "start solve" << endl;

	mlbm->Solve(20000, 200, 1);

	// for (int t=0; t<10000; ++t)
	// {
	// 	if (t%100==0)
	// 	{
	// 		cout << "Time= " << t << endl;
	// 		lbm->WriteFileH5(t,1);
	// 	}
	// 	lbm->CollideSRT();
	// 	lbm->ApplyWall();
	// 	lbm->Stream();
	// 	lbm->CalRhoV();	
	// }

	return 0;
}