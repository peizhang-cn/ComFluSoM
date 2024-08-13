/****************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics                *
 * Copyright (C) 2024 Pei Zhang                                             *
 * Email: peizhang.hhu@gmail.com                                            *
 *                                                                          *
 * This program is free software: you can redistribute it and/or modify     *
 * it under the terms of the GNU Affero General Public License as           *
 * published by the Free Software Foundation, either version 3 of the       *
 * License, or (at your option) any later version.                          *
 *                                                                          *
 * This program is distributed in the hope that it will be useful,          *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU Affero General Public License for more details.                      *
 *                                                                          *
 * You should have received a copy of the GNU Affero General Public License *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.   *
 * In cases where the constraints of the Open Source license prevent you 	*
 * from using ComFluSoM, please contact by peizhang.hhu@gmail.com for a 	*
 * commercial license. 														*
 ****************************************************************************/

#pragma once

// Functions for Metaball DEM

// read OBJ file data: vertices and faces, only works for triangle mesh now
// void ReadOBJ(string fname, vector<Vector3d>& ver, vector<VectorXi>& face)
// {
// 	string line;
// 	char c;
// 	int i, j, k;
// 	double x, y, z;
// 	string si, sj, sk;

// 	ifstream in(fname);
// 	if (!in)
// 	{
// 		std::cout << "Could not open file: " << fname << std::endl;
// 		abort();
// 	}
// 	while ( getline( in, line ) )                          
// 	{
// 		istringstream ss( line );
// 	  	if (line[0]=='v' && line[1]==' ')
// 	  	{
// 	  		ss >> c >> x >> y >> z;
// 	        Vector3d veri;
// 	        veri << x,y,z;
// 	        ver.push_back(veri);             
// 	  	}
// 	  	else if (line[0]=='f')
// 	  	{
// 			ss >> c >> si >> sj >> sk;                          
// 	        i = stoi(si)-1;  j = stoi(sj)-1;  k = stoi( sk )-1;
// 	        VectorXi facei;
// 	        facei.resize(3);
// 	        facei << i,j,k;
// 	        face.push_back(facei);  	
// 	  	}
// 	}
// 	in.close();
// }

void WriteOBJ(string fname, vector<Vector3d>& ver, vector<VectorXi>& face)
{
	ofstream outfile;
	outfile.open(fname);
	for (size_t i=0; i<ver.size(); ++i)
	{
		outfile << "v " << ver[i](0) << " " << ver[i](1) << " " << ver[i](2) << endl;
	}
	for (size_t i=0; i<face.size(); ++i)
	{
		outfile << "f " << face[i](0)+1 << "/1/1 " << face[i](1)+1<< "/1/1 " << face[i](2)+1<< "/1/1" << endl;
	}
}

// calculate meta function value
double CalMetaC(vector<Vector3d>& P, VectorXd& K, Vector3d& X)
{
	Vector3d W (1.,1.,1.);
	double c = 0.;
	for (size_t i=0; i<P.size(); ++i)
	{
		Vector3d xpi = X-P[i];
		Vector3d sqr = xpi.array()*xpi.array();
		double ri2 = W.dot(sqr);
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

// Find cross point between metaball and a vector, make sure that x0 is out of the metaball
double FindMetaVectorCrossPoints(vector<Vector3d>& P, VectorXd& K, double C, Vector3d& x0, Vector3d& vec, size_t inum, double tol)
{
	// x0: starting point of the vector, vec: direction, inum: maximum iteration number, tol: tolerance
	double c0;			// function value
	Vector3d f0;		// gradient
	Vector3d xi = x0;	// position
	double k = 0.;
	for (size_t i=0; i<inum; ++i)
	{
		CalMetaCF(P,K,xi,c0,f0);
		if (abs(C-c0)<tol)	break;

		k += (C-c0)/f0.dot(vec);	// first order taylor expension
		xi = x0+k*vec;
	}
	return k;
}

double FindMetaSegmentCrossPoints(vector<Vector3d>& P, VectorXd& K, double C, Vector3d& x0, Vector3d& x1)
{
	double k = 0.;
	double c0 = CalMetaC(P,K,x0);
	double c1 = CalMetaC(P,K,x1);
	k = (c0-C)/(c0-c1);
	return k;
}

double FindMetaUnitCubeCrossVolume(vector<Vector3d>& P, VectorXd& K, double C, Vector3d& x0)
{
	double c0;
	Vector3d f0;
	CalMetaCF(P,K,x0, c0, f0);
	f0.normalize();
	Vector3d xmin = x0 + 0.5*f0;
	Vector3d xmax = x0 - 0.5*f0;
	double cmin = CalMetaC(P,K,xmin);
	double cmax = CalMetaC(P,K,xmax);

	double k;
	if (C>cmax)			k = 0.;
	else if (C<cmin)	k = 1.;
	else				k = (cmax-C)/(cmax-cmin);
	return k;
}

double FindMetaUnitCubeCrossVolume(vector<Vector3d>& P, VectorXd& K, double C, double Cmin, double Cmax, Vector3d& x0, int n)
{
	double c0;
	Vector3d f0;
	CalMetaCF(P,K,x0, c0, f0);

	double v;
	if (c0<Cmin)		v = 0.;
	else if (c0>Cmax)	v = 1.;
	else
	{
	    int nc = 0.;
	    Vector3d xc (0.,0.,0.);

	    for (int ni=0; ni<n; ++ni)
	    for (int nj=0; nj<n; ++nj)
	    for (int nk=0; nk<n; ++nk)
	    {
	        xc(0) = x0(0)-0.5+(ni+0.5)/n;
	        xc(1) = x0(1)-0.5+(nj+0.5)/n;
	        xc(2) = x0(2)-0.5+(nk+0.5)/n;

	        double cc = CalMetaC(P,K,xc);
	        if (cc>C)    nc++;
	    }
		v = ((double) nc)/((double) n*n*n);
	}
	return v;
}

double FindMetaUnitCubeCrossVolume(vector<Vector3d>& P, VectorXd& K, double C, double Cmin, double Cmax, Vector3d& x0)
{
	double c0;
	Vector3d f0;
	CalMetaCF(P,K,x0, c0, f0);

	double v;
	if (c0<Cmin)		v = 0.;
	else if (c0>Cmax)	v = 1.;
	else
	{
		double dis = (C-c0)/f0.norm();
		if (dis<-0.5)		v = 1.;
		else if (dis<0.5)	v = 0.5-dis;
		else				v = 0.;	
	}
	return v;
}

// Make sure that x0 is within metaball
void FindMetaVectorCrossPointsWithRay(vector<Vector3d>& P, VectorXd& K, Vector3d& x0, Vector3d& vec, Vector3d& cp, size_t inum, double tol)
{
	// x0: starting point of the vector, vec: direction, sp: surface point
	Vector3d xi = x0;
	Vector3d xj = x0+vec;

	double c0 = CalMetaC(P, K, xi);
	double c1 = CalMetaC(P, K, xj);
	if (c0<1.)
	{
		cout << "x0 is out of the metaball!" << endl;
		abort();
	}
	// move xi and xj to make sure that ci>1 and cj<1.
	for (size_t s=1; s<1e5; ++s)
	{
		if (c1<1.)	break;
		else
		{
			xi = xj;
			xj = x0 + pow(2,s)*vec;
			c1 = CalMetaC(P, K, xj);
		}
	}
	// use middle point method to find crossing point
	for (size_t s=0; s<inum; ++s)
	{
		Vector3d xm = 0.5*(xi+xj);
		double cm = CalMetaC(P, K, xm);
		if (abs(cm-1.)<tol)
		{
			cp = xm;
			break;
		}
		if (cm>1.)	xi = xm;
		else		xj = xm;
	}
}


// get meta surface mesh by finding crossing points between every ray and meta, the ray is from unitSphere mesh
void GetMetaSurfaceMesh(vector<Vector3d>& unitSphere, vector<Vector3d>& P0, VectorXd K0, vector<Vector3d>& Pm, double rs)
{
	Pm.resize(unitSphere.size());
	// move all control points to mass center of control points (which should be inside of the meta for convex case)
	Vector3d mc0 = Vector3d::Zero(); 
	for (size_t i=0; i<P0.size(); ++i)	{mc0 += P0[i];}
	mc0 /= P0.size();

	// find crossing points
	for (size_t i=0; i< unitSphere.size(); ++i) 
	{
		FindMetaVectorCrossPointsWithRay(P0, K0, mc0, unitSphere[i], Pm[i], 1e5, 1.e-12);
		Pm[i] += rs*unitSphere[i];
	}
}

void BalanceSurfaceMesh(vector<VectorXi>& F, vector<Vector3d>& P, vector<Vector3d>& P0, VectorXd K0)
{
	double kn = 2.e-4;
	// double gn = 1.e-2;
	double dt = 0.1;
	double area = CalSurfaceArea(F, P)/F.size();
	double le = sqrt(4.*area/sqrt(3));

	size_t np = P.size();

	vector<Vector3d> forces(np);
	vector<Vector3d> vels(np);
	for (size_t i=0; i<np; ++i)
	{
		forces[i].setZero();
		vels[i].setZero();
	}

	double maxl=0.;
	double minl=1.e20;
	for (size_t i=0; i<F.size(); ++i)
	{
		size_t numf = F[i].size();
		for (size_t j=0; j<numf; ++j)
		{
			size_t ind0 = j;
			size_t ind1 = (j+1)%numf;
			Vector3d v0 = P[F[i](ind0)];				// the first vertice
			Vector3d v1 = P[F[i](ind1)];				// the second vertice

			double l = (v0-v1).norm();
			if (l<minl)	minl = l;
			if (l>maxl)	maxl = l;
		}		
	}

	cout << "maxl: " << maxl << endl;
	cout << "minl: " << minl << endl;

	for (size_t t=0; t<10000; ++t)
	{
		for (size_t i=0; i<F.size(); ++i)
		{
			size_t numf = F[i].size();
			for (size_t j=0; j<numf; ++j)
			{
				size_t ind0 = j;
				size_t ind1 = (j+1)%numf;
				Vector3d v0 = P[F[i](ind0)];				// the first vertice
				Vector3d v1 = P[F[i](ind1)];				// the second vertice

				double l = (v0-v1).norm();
				Vector3d n = (v0-v1)/l;
				Vector3d f = -kn*(l-le)/le*n;

				forces[F[i](ind0)] += f/* - gn*vels[F[i](ind0)]*/;
				forces[F[i](ind1)] -= f/* - gn*vels[F[i](ind1)]*/;
			}
		}

		for (size_t i=0; i<np; ++i)
		{
			Vector3d pc = P[i];
			double c0;
			Vector3d f0;
			CalMetaCF(P0,K0,pc,c0,f0);
			Vector3d n0 = f0;
			n0.normalize();
			vels[i] += forces[i]*dt;
			vels[i] -= vels[i].dot(n0)*f0;
			P[i] += vels[i]*dt;

			CalMetaCF(P0,K0,P[i],c0,f0);

			// cout << "c: " << c0 << endl;
			// f0.normalize();
			// P[i] += (1.-c0)*f0;
			P[i] += (1.-c0)/f0.squaredNorm()*f0;

			// if (i==0)
			// {
			// 	cout << "c: " << c0 << endl;
			// 	cout << "f0: " << f0.norm() << endl;
			// 	cout << "(1.-c0) " << (1.-c0) << endl;
			// }
			// cout << "f0: " << f0.norm() << endl;
			// cout << "(1.-c0) " << (1.-c0) << endl;
			// cout << "(1.-c0)*f0 " << (1.-c0)*f0.norm() << endl;
			CalMetaCF(P0,K0,P[i],c0,f0);
			// cout << "c: " << c0 << endl;

			// abort();
			// double tol = 1.e-8;
			// for (size_t m=0; m<1; ++m)
			// {
			// 	CalMetaCF(P0,K0,P[i],c0,f0);
			// 	f0.normalize();
			// 	P[i] += (1.-c0)*f0;
			// 	if (abs(1.-c0)<tol)	break;
			// }
			forces[i].setZero();
			vels[i] = P[i]-pc;
			// vels[i].setZero();

			// if (i==0)
			// {
			// 	cout << "after c: " << c0 << endl;
			// 	cout << "===========" << endl;
			// }
		}
	}

	maxl=0.;
	minl=1.e20;
	for (size_t i=0; i<F.size(); ++i)
	{
		size_t numf = F[i].size();
		for (size_t j=0; j<numf; ++j)
		{
			size_t ind0 = j;
			size_t ind1 = (j+1)%numf;
			Vector3d v0 = P[F[i](ind0)];				// the first vertice
			Vector3d v1 = P[F[i](ind1)];				// the second vertice

			double l = (v0-v1).norm();
			if (l<minl)	minl = l;
			if (l>maxl)	maxl = l;
		}
	}

	cout << "after" << endl;
	cout << "maxl: " << maxl << endl;
	cout << "minl: " << minl << endl;
	// abort();
}

void FindMetaInitPoint(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, Vector3d& initP)
{
	// find the closest control point pair to get the inital point for Newtons method.
	size_t p = 0;
	size_t q = 0;
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
}
// Loop over surface point list (Sp) and find the point that with maximum metaball function value (metaball defined by P0 and K0)
void FindClosestPointFromSurfacePoints(vector<Vector3d>& P0, VectorXd& K0, vector<Vector3d>& Sp, Vector3d& Cp)
{
	size_t id = 0;
	double cmax = 0.;
	for (size_t i=0; i<Sp.size(); ++i)
	{
		double c0 = CalMetaC(P0, K0, Sp[i]);
		if (c0>cmax)
		{
			cmax = c0;
			id = i;
		}
	}
	Cp = Sp[id];
}
// Use middle point of closest surface points as inital point
void FindMetaInitPointWithSurfacePoints(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, vector<Vector3d>& Sp0, vector<Vector3d>& Sp1, Vector3d& initP)
{
	Vector3d cp0, cp1;
	FindClosestPointFromSurfacePoints(P1, K1, Sp0, cp0);
	FindClosestPointFromSurfacePoints(P0, K0, Sp1, cp1);
	initP = 0.5*(cp0+cp1);
}

void FindMetaClosestPointsWithCell(size_t D, vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, size_t Level, double lsize, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1)
{
	vector<Vector3d> Ne;
	if (D==3)
	{
		Ne  = {     { 0, 0, 0},
					{ 0, 1, 0}, { 0, 0, 1}, {0, -1, 0}, { 0, 0,-1}, 
					{ 0, 1, 1}, {0, -1, 1}, {0, -1,-1}, { 0, 1,-1} };
	}
	else if (D==2)
	{
		Ne  = {		{ 0, 0, 0},
					{ 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
					{ 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0} };		
	}

	Vector3d xinit0;
	FindMetaInitPoint(P0, P1, K0, K1, xinit0);
	Vector3d xinit = xinit0;
	double cosa = 1.;
	double c0 = 0.;
	double c1 = 0.;
	Vector3d f0 (0.,0.,0.);
	Vector3d f1 (0.,0.,0.);
	for (size_t l=0; l<Level; ++l)
	{
		lsize /= pow(10,l+1);
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

	Xc0 = pk+(1.-c0)*f0;
	Xc1 = pk+(1.-c1)*f1;

	finalP = pk;
}

void FindMetaClosestPoints(vector<Vector3d>& P0, vector<Vector3d>& P1, VectorXd& K0, VectorXd& K1, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, bool& conv, size_t& cvt)
{
	Vector3d pk = initP;
	Vector3d dk = Vector3d::Zero();
	Vector3d dkn;

	Vector3d fp = Vector3d::Zero();
	Vector3d fq = Vector3d::Zero();
	// Newton's method
	for (size_t k=0; k<=100; ++k)
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

		PartialPivLU<MatrixXd> lu(j0+j1);
		dkn = -lu.inverse()*(f0+f1);

		double dif = abs(dkn.norm()-dk.norm());

		if (dif>1.e-12)
		{
			dk = dkn;
			pk += dkn;
		}
		else
		{
			fp = f0;
			fq = f1;
			cvt = k+1;
			break;		
		}
	}

	double c0 = CalMetaC(P0, K0, pk);
	double c1 = CalMetaC(P1, K1, pk);

	Xc0 = pk+(1.-c0)/fp.squaredNorm()*fp;
	Xc1 = pk+(1.-c1)/fq.squaredNorm()*fq;

	finalP = pk;

	if (cvt>90)	conv = false;
}

void FindMetaPlaneClosestPoints(Vector3d& N, Vector3d& X, vector<Vector3d>& P, VectorXd& K0, Vector3d& X0, bool first, Vector3d& initP, Vector3d& finalP, Vector3d& Xc0, Vector3d& Xc1, bool& conv, size_t& cvt)
{
	Vector3d ex (1.,0.,0.);								// x-axis
	Quaterniond q0 = Quaterniond::FromTwoVectors(N,ex); // rotation from N to x-axis
	Vector3d Xr = q0._transformVector(X-X0);			// rotate the point on plane
	vector<Vector3d> P0(P.size());						// points after rotation
	size_t ind = 0;											// the index of small distance control point
	double minDis = 1.0e30;
	for (size_t i=0; i<P.size(); ++i)
	{
		P0[i] = q0._transformVector(P[i]-X0);			// rotate all control points to current coordinate
		if (P0[i](0)<minDis)
		{
			ind = i;
			minDis = P0[i](0);
		}
	}
	Vector2d pk (initP(1), initP(2));
	if (first) pk << P0[ind](1), P0[ind](2);

	Vector2d dk = Vector2d::Zero();
	Vector2d dkn;	
	Vector2d fp;
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

		PartialPivLU<MatrixXd> lu(j0);
		dkn = -lu.inverse()*f0;

		if (abs(f0(0))<1.e-12 && abs(f0(1))<1.e-12)
		{
			fp = f0;
			cvt = k+1;
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

	double c;
	Vector3d ff;
	CalMetaCF(P0, K0, xc0r, c, ff);
	if (c<1.e-8)	conv = false;

	Vector3d xc1r = xc0r+(1.-c)/ff.squaredNorm()*ff;

	Quaterniond qi = q0.inverse();						// inverse rotation
	Xc0 = X0+qi._transformVector(xc0r);					// closest point on plane after rotated back
	Xc1 = X0+qi._transformVector(xc1r);					// closest point on metaball after rotated back
}