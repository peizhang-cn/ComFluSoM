// Functions for Metaball DEM

// read OBJ file data: vertices and faces, only works for triangle mesh now
void ReadOBJ(string fname, vector<Vector3d>& ver, vector<VectorXi>& face)
{
	string line;
	char c;
	int i, j, k;
	double x, y, z;
	string si, sj, sk;

	ifstream in(fname);
	if (!in)
	{
		std::cout << "Could not open file: " << fname << std::endl;
		abort();
	}
	while ( getline( in, line ) )                          
	{
		istringstream ss( line );
	  	if (line[0]=='v' && line[1]==' ')
	  	{
	  		ss >> c >> x >> y >> z;
	        Vector3d veri;
	        veri << x,y,z;
	        ver.push_back(veri);             
	  	}
	  	else if (line[0]=='f')
	  	{
			ss >> c >> si >> sj >> sk;                          
	        i = stoi(si)-1;  j = stoi(sj)-1;  k = stoi( sk )-1;
	        VectorXi facei;
	        facei.resize(3);
	        facei << i,j,k;
	        face.push_back(facei);  	
	  	}
	}
	in.close();
}

void WriteOBJ(string fname, vector<Vector3d>& ver, vector<VectorXi>& face)
{
	ofstream outfile;
	outfile.open(fname);
	for (size_t i=0; i<ver.size(); ++i)
	{
		outfile << "v " << ver[i](0) << " " << ver[i](1) << " " << ver[i](2) << endl;
	}
	// for (size_t i=0; i<face.size(); ++i)
	// {
	// 	size_t a = face[i](0);
	// 	size_t b = face[i](1);
	// 	size_t c = face[i](2);
	// 	Vector3d n = (ver[b]-ver[a]).cross(ver[c]-ver[a]);
	// 	n.normalize();
	// 	outfile << "vn " << n(0) << " " << n(1) << " " << n(2) << endl;
	// }
	for (size_t i=0; i<face.size(); ++i)
	{
		// outfile << "f " << face[i](0)+1 << "/"<<i <<"/1" << " " << face[i](1)+1<< "/"<<i<<"/1" << " " << face[i](2)+1<< "/"<<i<<"/1" << endl;
		outfile << "f " << face[i](0)+1 << "/1/1 " << face[i](1)+1<< "/1/1 " << face[i](2)+1<< "/1/1" << endl;
	}
}

// calculate meta function value
double CalMetaC(vector<Vector3d>& P, VectorXd& K, Vector3d& X)
{
	double c = 0.;
	for (size_t i=0; i<P.size(); ++i)
	{
		double ri2 = (X-P[i]).squaredNorm();
		c += K(i)/ri2;
	}
	return c;
}

void CalMetaCF(vector<Vector3d>& P, VectorXd& K, Vector3d& X, double& C, Vector3d& F)
{
	double c = 0.;
	Vector3d f = Vector3d::Zero();
	for (size_t i=0; i<P.size(); ++i)
	{
		Vector3d a = X-P[i];
		double ri2inv = 1./a.squaredNorm();
		double ri4inv = ri2inv*ri2inv;
		f -= 2.*K(i)*ri4inv*a;
		c += K(i)*ri2inv;
	}
	C = c;
	F = f;
}


// Assuming that the origin is inside of meta and close to mass center, find the crossing point between meta and the ray (from origin with dirction N)
double CalMetaRayCrossPoint(Vector3d& N, vector<Vector3d>& P0, VectorXd& K0, double r0)
{
	double n2 = N.squaredNorm();
	n2 = 1.;
	vector<double> p2, pn;
	p2.resize(P0.size());
	pn.resize(P0.size());
	for (size_t i=0; i<P0.size(); ++i)
	{
		p2[i] = P0[i].squaredNorm();
		pn[i] = -2.*P0[i].dot(N);
	}
	// Newton's method
	double r = r0;
	for (size_t k=0; k<10000; ++k)
	{
		double f = -1.;
		double df = 0.;
		for (size_t i=0; i<P0.size(); ++i)
		{
			double b = (n2*r*r+pn[i]*r+p2[i]);
			f += K0(i)/b;
			df += -K0(i)/(b*b)*(2.*n2*r+pn[i]);
		}
		r -= f/df;
		if (f<1.0e-12)	break;
	}
	return r;
}

Vector3d CalMetaRayCrossPoint(Vector3d& N, vector<Vector3d>& P0, VectorXd& K0)
{
	Vector3d x = N;
	for (size_t i=2; i<1000; ++i)
	{
		double f = CalMetaC(P0, K0, x);
		if (f<1.)	break;
		else		x *= i;
	}
	Vector3d x0 = Vector3d::Zero();
	Vector3d x1 = x;
	Vector3d xm;
	for (size_t i=0; i<1000; ++i)
	{
		xm = 0.5*(x0+x1);
		double f = CalMetaC(P0, K0, xm);
		if (abs(f-1.)<1.0e-12)	break;
		else if (f<1.)			x1 = xm;
		else					x0 = xm;
	}
	return xm;
}

// get meta surface mesh by finding crossing points between every ray and meta, the ray is from unitSphere mesh
void GetMetaSurfaceMesh(vector<Vector3d>& unitSphere, vector<Vector3d>& P0, VectorXd K0, vector<Vector3d>& Pm)
{
	Pm.resize(unitSphere.size());
	// move all control points to mass center of control points (which should be inside of the meta for convex case)
	Vector3d mc0 = Vector3d::Zero(); 
	for (size_t i=0; i<P0.size(); ++i)	{mc0 += P0[i];}
	mc0 /= P0.size();
	// control points' relative position to mc0
	vector<Vector3d> P;
	P.resize(P0.size());
	for (size_t i=0; i<P0.size(); ++i)	{P[i] = P0[i]-mc0;}
	// find crossing points
	for (size_t i=0; i< unitSphere.size(); ++i) 
	{
		double r = CalMetaRayCrossPoint(unitSphere[i], P, K0, 1.);
		if (r<0.)
		{
			cout << "r in : " << r << endl;
			cout << "negative!" << endl;
			abort();
		}
		Pm[i] = mc0 + r*unitSphere[i];

		Pm[i] = mc0 + CalMetaRayCrossPoint(unitSphere[i], P, K0);
	}
}

void FindMetaInitPoint(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, Vector3d& initP)
{
	// find the closest control point pair to get the inital point for Newtons method.
	size_t p,q;
	double minDis2 = 1.0e30;
	for (size_t i=0; i<P0.size(); ++i)
	for (size_t j=i; j<P1.size(); ++j)
	{
		double dis2 = (P0[i]-P1[j]).squaredNorm();
		if (dis2<minDis2)
		{
			p = i;
			q = j;
			minDis2 = dis2;
		}
	}
	// Inital point
	Vector3d xinit0 = 0.5*(P0[p]+P1[q]);

	Vector3d p0 = P0[p];
	Vector3d p1 = P1[q];

	for (size_t i=0; i<10000; ++i)
	{
		double cm0 = CalMetaC(P0, K0, xinit0);
		double cm1 = CalMetaC(P1, K1, xinit0);
		if (cm0<1. && cm1<1.)	break;
		else if (cm0>1.)		p0 = xinit0;
		else if (cm1>1.)		p1 = xinit0;
		xinit0 = 0.5*(p0+p1);
	}
	initP = xinit0;
	// cout << "xinit: " << xinit.transpose() << endl;
	// cout << "xinit: " << xinit.transpose() << endl;
}

void FindMetaClosestPointsWithCell(vector<Vector3d>& Ne, vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, size_t Level, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1)
{
	Vector3d xinit0;
	// FindMetaInitPoint(Ne, P0, P1, K0, K1, Level, pk);
	FindMetaInitPoint(P0, P1, K0, K1, xinit0);

	Vector3d xinit = xinit0;
	double lsize = 0.001;
	double cosa = 1.;
	double c0,c1;
	Vector3d f0, f1;
	for (size_t l=0; l<Level; ++l)
	{
		lsize /= pow(10,l);
		for (size_t s=0; s<1000000; ++s)
		{
			Vector3d xinits = xinit;
			for (size_t n=0; n<Ne.size(); ++n)
			{
				Vector3d xn = xinit + lsize*Ne[n];
				CalMetaCF(P0, K0, xn, c0, f0);
				CalMetaCF(P1, K1, xn, c1, f1);
				double cosn = f0.dot(f1)/(f0.norm()*f1.norm());
				if ((cosn < cosa) && (c0<1.) && (c1<1.))
				{
					xinit = xn;
					cosa = cosn;
				}
			}
			if (xinits == xinit)	break;
		}
		if (abs(1.+cosa)<1.e-12) 	break;
	}
	Vector3d pk = xinit;

	// cout << "cosa: " << cosa << endl;
 
	// // double c0, c1;
	// Vector3d fp, fq;
	// CalMetaCF(P0, K0, pk, c0, fp);
	// CalMetaCF(P1, K1, pk, c1, fq);
	// 	cout << "c0: " << c0 << endl;
	// 	cout << "c1: " << c1 << endl;
	// 	cout << "fp: " << fp.transpose() << endl;
	// 	cout << "fq: " << fq.transpose() << endl;
	// 	cout << "++++++++++" << endl;
	// // contact = false;
	Xc0 = pk+(1.-c0)*f0;
	Xc1 = pk+(1.-c1)*f1;

	// double c00 = CalMetaC(P0, K0, Xc0);
	// double c11 = CalMetaC(P1, K1, Xc1);

	// Xc0 = pk+sqrt(1./CalMetaC(P0, K0, Xc0))*(1.-c0)*f0;
	// Xc1 = pk+sqrt(1./CalMetaC(P1, K1, Xc1))*(1.-c1)*f1;
	finalP = pk;
}

void FindMetaPlaneClosestPointsWithCell(Vector3d& N, Vector3d& X, vector<Vector3d>& P, VectorXd& K, Vector3d& X0, size_t Level, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1)
{
	Vector3d ex (1.,0.,0.);								// x-axis
	Quaterniond q0 = Quaterniond::FromTwoVectors(N,ex); // rotation from N to x-axis
	Vector3d Xr = q0._transformVector(X-X0);			// rotate the point on plane
	vector<Vector3d> P0(P.size());						// points after rotation
	size_t ind;											// the index of small distance control point
	double minDis = 1.0e30;
	for (size_t i=0; i<P.size(); ++i)
	{
		P0[i] = q0._transformVector(P[i]-X0);
		if (P0[i](0)<minDis)
		{
			ind = i;
			minDis = P0[i](0);
		}
	}

	vector<Vector3d> Ne;
	Ne  = {     { 0, 0, 0},
				{ 0, 1, 0}, { 0, 0, 1}, {0, -1, 0}, { 0, 0,-1}, 
				{ 0, 1, 1}, {0, -1, 1}, {0, -1,-1}, { 0, 1,-1} };

	// Inital point
	Vector3d xinit0 (Xr(0), P0[ind](1), P0[ind](2));
	Vector3d xinit = xinit0;
	double lsize = 0.01;
	double cosa = 1.;
	double c0;
	Vector3d f0;
	for (size_t l=0; l<Level; ++l)
	{
		lsize /= pow(10,l);
		for (size_t s=0; s<1000000; ++s)
		{
			Vector3d xinits = xinit;
			for (size_t n=0; n<Ne.size(); ++n)
			{
				Vector3d xn = xinit + lsize*Ne[n];
				CalMetaCF(P0, K, xn, c0, f0);
				double cosn = -f0(0)/f0.norm();
				if ((cosn < cosa) && (c0<1.))
				{
					xinit = xn;
					cosa = cosn;
				}
			}
			if (xinits == xinit)	break;
		}
	}
	// if (abs(1.+cosa)>1.0e-8)
	// {
	// 	cout << "cosa: " << cosa << endl;
	// 	abort();
	// }
	finalP = xinit;
	Vector3d xc0r = finalP;
	Vector3d xc1r = xc0r;
	xc1r(0) += (1.-c0);

	Quaterniond qi = q0.inverse();						// inverse rotation
	Xc0 = X0+qi._transformVector(xc0r);					// closest point on plane after rotated back
	Xc1 = X0+qi._transformVector(xc1r);					// closest point on metaball after rotated back
}


void FindMetaClosestPoints(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, bool& conv)
{
	Vector3d pk = initP;
	Vector3d dk = Vector3d::Zero();
	Vector3d dkn;

	Vector3d fp = Vector3d::Zero();
	Vector3d fq = Vector3d::Zero();
	// Matrix3d jp, jq;
	// Newton's method
	size_t cvt = 0;
	for (size_t k=0; k<=10000; ++k)
	{
		Vector3d f0 = Vector3d::Zero();
		Vector3d f1 = Vector3d::Zero();
		Matrix3d j0 = Matrix3d::Zero();
		Matrix3d j1 = Matrix3d::Zero();
		for (size_t i=0; i<P0.size(); ++i)
		{
			Vector3d a = pk-P0[i];
			double ri2inv = 1./a.squaredNorm();
			double ri4inv = ri2inv*ri2inv;
			double ri6inv = ri4inv*ri2inv;
			f0 -= 2.*K0(i)*ri4inv*a;
			j0(0,0) += 2.*K0(i)*(-ri4inv + 4.*a(0)*a(0)*ri6inv);
			j0(1,1) += 2.*K0(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
			j0(2,2) += 2.*K0(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
			j0(0,1) += 8.*K0(i)*a(0)*a(1)*ri6inv;
			j0(0,2) += 8.*K0(i)*a(0)*a(2)*ri6inv;
			j0(1,2) += 8.*K0(i)*a(1)*a(2)*ri6inv;
		}
		j0(1,0) = j0(0,1);
		j0(2,0) = j0(0,2);
		j0(2,1) = j0(1,2);
		for (size_t i=0; i<P1.size(); ++i)
		{
			Vector3d a = pk-P1[i];
			double ri2inv = 1./a.squaredNorm();
			double ri4inv = ri2inv*ri2inv;
			double ri6inv = ri4inv*ri2inv;
			f1 -= 2.*K1(i)*ri4inv*a;
			j1(0,0) += 2.*K1(i)*(-ri4inv + 4.*a(0)*a(0)*ri6inv);
			j1(1,1) += 2.*K1(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
			j1(2,2) += 2.*K1(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
			j1(0,1) += 8.*K1(i)*a(0)*a(1)*ri6inv;
			j1(0,2) += 8.*K1(i)*a(0)*a(2)*ri6inv;
			j1(1,2) += 8.*K1(i)*a(1)*a(2)*ri6inv;
		}
		j1(1,0) = j1(0,1);
		j1(2,0) = j1(0,2);
		j1(2,1) = j1(1,2);

		dkn = -(j0+j1).inverse()*(f0+f1);
		double dif = abs(dkn.norm()-dk.norm());

		if (dif>1.e-12)
		{
			dk = dkn;
			pk += dkn;
		}
		else
		{
			fp = f0.normalized();
			fq = f1.normalized();
			// fp = f0;
			// fq = f1;
			// jp = j0;
			// jq = j1;
			cvt = k;
			break;		
		}
	}
	// cout << "cvt: " << cvt << endl;
	// cout << "pk: " << pk.transpose() << endl;
	double c0 = CalMetaC(P0, K0, pk);
	double c1 = CalMetaC(P1, K1, pk);
	// 	cout << "c0: " << c0 << endl;
	// 	cout << "c1: " << c1 << endl;
	// // contact = false;
	Xc0 = pk+(1.-c0)*fp;
	Xc1 = pk+(1.-c1)*fq;

	// double c00 = CalMetaC(P0, K0, Xc0);
	// double c11 = CalMetaC(P1, K1, Xc1);

	// Xc0 = pk+sqrt(1./CalMetaC(P0, K0, Xc0))*(1.-c0)*fp;
	// Xc1 = pk+sqrt(1./CalMetaC(P1, K1, Xc1))*(1.-c1)*fq;
	finalP = pk;

	if (cvt==100)	conv = false;
	double cosn = fp.dot(fq);
	if (abs(1.+cosn)>1.e12)
	{
		cout << "abs(1.+cosn)>1.e12 : " << cosn << endl; 
		abort();
	}
	// if (finalP.norm()>1e8) conv = false;
	// if (finalP.norm()>1e8)
	// {
	// 	cout << "initP: " << initP.transpose() << endl;
	// 	cout << "finalP: " << finalP.transpose() << endl;
	// 	cout << (Xc1-Xc0).norm() << endl;
	// 	cout << "Xc0: " << Xc0.transpose() << endl;
	// 	cout << "Xc1: " << Xc1.transpose() << endl;
	// 	cout << "c0: " << c0 << endl;
	// 	cout << "c1: " << c1 << endl;
	// 	cout << CalMetaC(P0, K0, initP) << endl;
	// 	cout << CalMetaC(P1, K1, initP) << endl;
	// 	// abort();
	// }
}

void FindMetaPlaneClosestPoints(Vector3d& N, Vector3d& X, vector<Vector3d>& P, VectorXd& K0, Vector3d& X0, bool first, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, bool& conv)
{
	Vector3d ex (1.,0.,0.);								// x-axis
	Quaterniond q0 = Quaterniond::FromTwoVectors(N,ex); // rotation from N to x-axis
	Vector3d Xr = q0._transformVector(X-X0);			// rotate the point on plane
	vector<Vector3d> P0(P.size());						// points after rotation
	size_t ind;											// the index of small distance control point
	double minDis = 1.0e30;
	for (size_t i=0; i<P.size(); ++i)
	{
		P0[i] = q0._transformVector(P[i]-X0);
		if (P0[i](0)<minDis)
		{
			ind = i;
			minDis = P0[i](0);
		}
	}
	Vector2d pk (initP(1), initP(2));
	if (first) pk << P0[ind](1), P0[ind](2);

	Vector2d initP0 = pk;
	Vector2d dk = Vector2d::Zero();
	Vector2d dkn;	
	Vector2d fp;
	size_t cvt = 0;
	for (size_t k=0; k<10000; ++k)
	{
		Vector2d f0 = Vector2d::Zero();
		Matrix2d j0 = Matrix2d::Zero();
		for (size_t i=0; i<P0.size(); ++i)
		{
			Vector3d pki (Xr(0), pk(0), pk(1));
			Vector3d a = pki-P0[i];
			double ri2inv = 1./a.squaredNorm();
			double ri4inv = ri2inv*ri2inv;
			double ri6inv = ri4inv*ri2inv;
			f0(0) -= 2.*K0(i)*ri4inv*a(1);
			f0(1) -= 2.*K0(i)*ri4inv*a(2);

			j0(0,0) += 2.*K0(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
			j0(1,1) += 2.*K0(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
			j0(0,1) += 8.*K0(i)*a(1)*a(2)*ri6inv;
		}
		j0(1,0) = j0(0,1);

		dkn = -j0.inverse()*f0;

		if (abs(f0(0))<1.e-12 && abs(f0(1))<1.e-12)
		{
			fp = f0;
			cvt = k;
			break;		
		}
		else
		{
			dk = dkn;
			pk += dkn;
		}
	}
	finalP << 0., pk(0), pk(1);
	Vector3d xc0r (Xr(0), pk(0), pk(1));				// closest point on the plane after rotation
	// double c = CalMetaC(P0, K0, xc0r);
	double c;
	Vector3d ff;
	CalMetaCF(P0, K0, xc0r, c, ff);
	if (c<1.e-8)	conv = false;

	// if (ff(0)<0.)
	// {
	// 	conv = false;
	// 	cout << "ff: " << ff.transpose() << endl;
	// 	// abort();
	// 	// conv = false;
	// }

	Vector3d xc1r = xc0r;
	xc1r(0) += (1.-c);

	Quaterniond qi = q0.inverse();						// inverse rotation
	Xc0 = X0+qi._transformVector(xc0r);					// closest point on plane after rotated back
	Xc1 = X0+qi._transformVector(xc1r);					// closest point on metaball after rotated back

	double delta = 0.2-(Xc0-Xc1).norm();

	// if (delta>0.0002)
	// {
	// 	CalMetaCF(P0, K0, xc1r, c, ff);
	// 	cout << "c: " << c << endl;
	// 	cout << "ff: " << ff.transpose() << endl;
	// }

	// cout << "N: " << N.transpose() << endl;
	// cout << "Xc0: " << Xc0.transpose() << endl;
	// cout << "Xc1: " << Xc1.transpose() << endl;
	// cout << "Xc0r: " << xc0r.transpose() << endl;
	// cout << "Xc1r: " << xc1r.transpose() << endl;
	// cout << "delta: " << (Xc0-Xc1).norm() << endl;
	// cout << "cvt: " << cvt << endl;
	// abort();
	// if (cvt>100)
	// {
	// 	cout << "cvt>100" << endl;
	// 	dk.setZero();
	// 	pk = initP0;
	// 	if (first)	cout << "first!!!!" << endl;
	// 	cout << "initP0: " << initP0.transpose() << endl;
	// 	cout << "P0[ind]: " << P0[ind].transpose() << endl;

	// 	for (size_t k=0; k<10000; ++k)
	// 	{
	// 		Vector2d f0 = Vector2d::Zero();
	// 		Matrix2d j0 = Matrix2d::Zero();
	// 		for (size_t i=0; i<P0.size(); ++i)
	// 		{
	// 			Vector3d pki (Xr(0), pk(0), pk(1));
	// 			Vector3d a = pki-P0[i];
	// 			double ri2inv = 1./a.squaredNorm();
	// 			double ri4inv = ri2inv*ri2inv;
	// 			double ri6inv = ri4inv*ri2inv;
	// 			f0(0) -= 2.*K0(i)*ri4inv*a(1);
	// 			f0(1) -= 2.*K0(i)*ri4inv*a(2);

	// 			j0(0,0) += 2.*K0(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
	// 			j0(1,1) += 2.*K0(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
	// 			j0(0,1) += 8.*K0(i)*a(1)*a(2)*ri6inv;
	// 		}
	// 		j0(1,0) = j0(0,1);

	// 		dkn = -j0.inverse()*f0;

	// 		cout << "f: " << f0.transpose() << endl;
	// 		cout << "j: " << endl;
	// 		cout << j0 << endl;

	// 		if (abs(f0(0))<1.e-12 && abs(f0(1))<1.e-12)
	// 		{
	// 			fp = f0;
	// 			cvt = k;
	// 			break;		
	// 		}
	// 		else
	// 		{
	// 			dk = dkn;
	// 			pk += dkn;
	// 			cout << "N: " << N.transpose() << endl;
	// 			cout << "pki: " << pk.transpose() << endl;
	// 		}
	// 	}
	// 	cout << "pk111: " << pk.transpose() << endl;
	// 	abort();
	// }
}

// NOT WORKING
// Only works for cylinder drum. The axis of drum is x-axis as defult
void FindMetaDrumClosestPoints(Vector3d& Xd, double Rd, vector<Vector3d>& P0, VectorXd& K0, bool first, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, bool& conv)
{
	Vector3d pk = initP;
	Vector3d dk = Vector3d::Zero();
	Vector3d dkn;

	Vector3d fp = Vector3d::Zero();
	Vector3d fq = Vector3d::Zero();

	if (first)
	{
		size_t ind;
		double maxDis = -1.;
		for (size_t i=0; i<P0.size(); ++i)
		{
			Vector3d vri = P0[i]-Xd;
			vri(0) = 0.;
			double ri = vri.norm();
			if (vri.norm()>maxDis)
			{
				ind = i;
				maxDis = ri;
			}
		}
		Vector3d vr = P0[ind]-Xd;	// relative vector from drum center to pk
		vr(0) = 0.;
		pk = P0[ind] + vr.norm()/Rd*vr;
	}
	// Newton's method
	size_t cvt = 0;
	for (size_t k=0; k<=10000; ++k)
	{
		Vector3d f0 = Vector3d::Zero();
		Vector3d f1 = Vector3d::Zero();
		Matrix3d j0 = Matrix3d::Zero();
		Matrix3d j1 = Matrix3d::Zero();
		// for metaball
		for (size_t i=0; i<P0.size(); ++i)
		{
			Vector3d a = pk-P0[i];
			double ri2inv = 1./a.squaredNorm();
			double ri4inv = ri2inv*ri2inv;
			double ri6inv = ri4inv*ri2inv;
			f0 -= 2.*K0(i)*ri4inv*a;
			j0(0,0) += 2.*K0(i)*(-ri4inv + 4.*a(0)*a(0)*ri6inv);
			j0(1,1) += 2.*K0(i)*(-ri4inv + 4.*a(1)*a(1)*ri6inv);
			j0(2,2) += 2.*K0(i)*(-ri4inv + 4.*a(2)*a(2)*ri6inv);
			j0(0,1) += 8.*K0(i)*a(0)*a(1)*ri6inv;
			j0(0,2) += 8.*K0(i)*a(0)*a(2)*ri6inv;
			j0(1,2) += 8.*K0(i)*a(1)*a(2)*ri6inv;
		}
		j0(1,0) = j0(0,1);
		j0(2,0) = j0(0,2);
		j0(2,1) = j0(1,2);
		// for drum
		double r2i = 2./(Rd*Rd);
		f1(1) = r2i*(pk(1)-Xd(1));
		f1(2) = r2i*(pk(2)-Xd(2));
		j1(1,1) = j1(2,2) = r2i;

		dkn = -(j0+j1).inverse()*(f0+f1);
		double dif = abs(dkn.norm()-dk.norm());

		if (dif>1.e-12)
		{
			dk = dkn;
			pk += dkn;
		}
		else
		{
			fp = f0.normalized();
			cvt = k;
			break;
		}

		Vector3d vr = pk-Xd;	// relative vector from drum center to pk
		vr(0) = 0.;
		Xc0 = pk + vr.norm()/Rd*vr;

		double c = CalMetaC(P0, K0, pk);
		Xc1 = pk+(1.-c)*fp;

		cout << "pk: " << pk.transpose() << endl;

		cout << "Xc0: " << Xc0.transpose() << endl;
		cout << "Xc1: " << Xc1.transpose() << endl;
		abort();
	}
}