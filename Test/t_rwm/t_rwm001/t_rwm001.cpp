#include <RWM.h>

int main(int argc, char const *argv[])
{
	int nx = 50;
	int ny = 50;
	int nz = 0;
	RWM* a = new RWM(nx, ny, nz, 0.01, 1.);
	a->Init();
	a->Nproc = 12;
	Vector3d pos (0.5*nx, 0.5*ny, 0.);
	for (size_t p=0; p<100000; ++p)
	{
		a->AddParticle(-1,pos, 1.);
	}
	// double r = 10;

	// d->Move(1.);

	// for (int p=0; p<d->Lp.size(); ++p)
	// {
	// 	// cout << p << endl;
	// 	// cout << d->Lp[p]->X(0) << endl;
	// 	double dis = GetDistance(d->Lp[p]->X);
	// 	if (dis>r) cout << "wrong dis= " << dis << endl;
	// }
	// cout << d->Lp[999999]->X(0) << endl;
	// cout << d->Lp[999999]->X(1) << endl;
	// cout << d->Lp[999999]->X(2) << endl;
/*
	Vector_Field v, v1;
	Vector_Field* v_ptr;

	Vec3_t initV (0,0,0);
	v.resize(Nx+1, vector< vector<Vec3_t> >(Ny+1,vector<Vec3_t>(Nz+1,initV)));
*/
	// v_ptr = &v;

	// for (int i = 0; i <= Nx; ++i)
 //    for (int j = 0; j <= Ny; ++j)
 //    for (int k = 0; k <= Nz; ++k)
 //    {
 //    	v[i][j][k]= i, j, k;
 //    }

/*    d->V_ptr = &v;*/

	// Vec3_t x0 (13.2, 14.6, 18.9);

	// Vec3_t v0;

 //    d->CalV(x0, d->V_ptr, v0);

	Vector3d xs (25, 35, 0);
	Vector3d xsb;
	double r = 5;
	Vector3d v (0.,-0.005, 0.);
	for (int t=0; t<100000; ++t)
	{
		xsb = xs;
		xs += v;
		// cout << "move 1" << endl;
		a->Move();
		for (size_t p=0; p<a->Lp.size(); ++p)
		{
			if ((a->Lp[p]->X-xs).norm()<r)
			{
				// cout << "start" << endl;
				// double A = (a->Lp[p]->X-a->Lp[p]->Xb).squaredNorm();
				// double B = 2*(a->Lp[p]->Xb-xs).dot(a->Lp[p]->X-a->Lp[p]->Xb);
				// double C = (a->Lp[p]->Xb-xs).squaredNorm()-r*r;
				// double delta = sqrt(B*B-4*A*C);
				// double k0 = (-B+delta)/(2*A);
				// double k1 = (-B-delta)/(2*A);
				// double k = k0;
				// if (k1>=0.&&k1<=1.)		k=k1;
				// Vector3d xc = a->Lp[p]->Xb + k*(a->Lp[p]->X-a->Lp[p]->Xb);
				// Vector3d n = (xc-xs).normalized();
				// cout << "===================" << endl;
				// cout << "A= " << A << endl;
				// cout << "B= " << B << endl;
				// cout << "C= " << C << endl;
				// cout << "delta= " << B*B-4*A*C << endl;
				// cout << "dis0= " << (a->Lp[p]->X-xs).norm() << endl;
				// cout << "dis1= " << (a->Lp[p]->Xb-xs).norm() << endl;
				// cout << "k0= "<< k0 << endl;
				// cout << "k1= "<< k1 << endl;
				// cout << "k= "<< k << endl;
				// cout << "x= " << a->Lp[p]->X.transpose() << endl;
				// cout << "xb= " << a->Lp[p]->Xb.transpose() << endl;
				// cout << "xc= " << xc.transpose() << endl;
				// cout << "n= " << n.transpose() << endl;
				Vector3d xc, n;
				a->ContactPointWithSphere(a->Lp[p], xs, xsb, r, xc, n);
				a->Reflection(a->Lp[p], n, xc);
				// cout << "reflect x= " << a->Lp[p]->X.transpose() << endl;
				// abort();
				if ((a->Lp[p]->X-xs).norm()<r)
				{
					a->Lp[p]->X = xs+(r+1.0e-12)*(a->Lp[p]->X-xs).normalized();
					// cout << "still inside" << endl;
					// cout << "A= " << A << endl;
					// cout << "B= " << B << endl;
					// cout << "C= " << C << endl;
					// cout << "delta= " << B*B-4*A*C << endl;
					// cout << "dis0= " << (a->Lp[p]->X-xs).norm() << endl;
					// cout << "dis1= " << (a->Lp[p]->Xb-xs).norm() << endl;
					// cout << "k0= "<< k0 << endl;
					// cout << "k1= "<< k1 << endl;
					// cout << "k= "<< k << endl;
					// cout << "x= " << a->Lp[p]->X.transpose() << endl;
					// cout << "xc= " << xc.transpose() << endl;
					// cout << "n= " << n.transpose() << endl;
					// abort();
				}
				if ((a->Lp[p]->X-xs).norm()<r)
				{
					cout << "still inside" << endl;
					cout << "dis0= " << (a->Lp[p]->X-xs).norm() << endl;
					abort();
				}
			}
		}

		a->RemoveParticles();
		// d->Boundary();
		// cout << "move 2" << endl;
		if (t%100==0)
		{
			cout << "time step= " << t << endl;
			a->CalC();
			a->WriteFileH5(t);
		}
	}
    
    // cout << d->C[50][50][0] << endl;
	// Vec3_t v0, vx;

	// v0 = 20.5665, 14.21322, 37.1984;

	// d->CalV(v0, d->V_ptr, vx);

	// d->CalV(d->Lp[10]->X, v_ptr, vx);

	// cout << d->Lp[10]->X << endl;
	// cout << vx << endl;
	// v_ptr = &v;

	// v1 = *v_ptr;

	// cout <<v[50][50][50] << endl;
	// cout <<(*v_ptr)[50][50][50] << endl;

	// vector <double> s(10,1.2);

	// vector <double>* s_ptr;

	// s_ptr = &s;

	// double ss = s_ptr->at(1);

	// cout << ss << endl;

	return 0;
}