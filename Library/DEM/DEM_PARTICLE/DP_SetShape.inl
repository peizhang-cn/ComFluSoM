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

#ifndef DEM_PARTICLE_SET_SHAPE_H
#define DEM_PARTICLE_SET_SHAPE_H

inline void DEM_PARTICLE::RerangeFaceElementOrder() // only works for convex shape
{
	// rerange face elements order (right hand role for outside normal)
	for (size_t i=0; i<Faces.size();++i)
	{
		size_t ind;	// index of the point that not belongs to this face i
		for (size_t j=0; j<P0.size(); ++j)
		{
			bool exist = false;
			for (int k=0; k<Faces[i].size();++k)
			{
				if (j-Faces[i](k)==0)	// check if j belongs to face i
				{
					exist = true;
					break;
				}
			}
			if (!exist)
			{
				ind = j;
				break;
			}
		}
		Vector3d v0 = P0[Faces[i](0)]-P0[Faces[i](1)];
		Vector3d v1 = P0[Faces[i](2)]-P0[Faces[i](1)];
		Vector3d v3 = P0[ind]-P0[Faces[i](1)];
		if ((v1.cross(v0)).dot(v3)>0)	// check norm is with same side of point ind
		{
			VectorXi facer = Faces[i].reverse();	// if so reverse order of face i
			Faces[i] = facer;
		}
	}
}

inline void DEM_PARTICLE::SetPolygon3D(vector<Vector3d> ver0)
{
	ShapeType = 0;
	X.setZero();
	M = 1.0e32;
	I(0) = I(1) = I(2) = 1.0e12;
	size_t nv = ver0.size();	// number of vertices
	// faces
	VectorXi face(nv);
	face << 0, 1, 2, 3;
	Faces.push_back(face);
	Nfe = 5;
	// assuming face is convex and anti clockwise order
	Vector3d norm = (ver0[Faces[0](1)] - ver0[Faces[0](0)]).cross(ver0[Faces[0](2)] - ver0[Faces[0](1)]);
	Normal = norm;
	Normal.normalize();

	Vector3d n0 (0.,0.,1.);
	Quaterniond q0;
	q0.setFromTwoVectors(n0,Normal);

	Q0 = q0;
	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();
	// rotate all vertices
	vector<Vector3d> ver(nv);
	for (size_t i=0; i<nv; ++i)
	{
		size_t j = (i+1)%nv;
		Vector2i edge (i, j);
		Edges.push_back(edge);
		ver[i] = Qfi._transformVector(ver0[i]);
	}

	double cx = 0.; double cy = 0.;
	double jx = 0.; double jy = 0.;
	double sa = 0.;
	for (size_t i=0; i<nv; ++i)
	{
		size_t j = (i+1)%nv;
		double fact = (ver[i](0)*ver[j](1) - ver[j](0)*ver[i](1));
		cx += (ver[i](0)+ver[j](0))*fact;
		cy += (ver[i](1)+ver[j](1))*fact;
		jx += (ver[i](1)*ver[i](1) + ver[i](1)*ver[j](1) + ver[j](1)*ver[j](1))*fact;
		jy += (ver[i](0)*ver[i](0) + ver[i](0)*ver[j](0) + ver[j](0)*ver[j](0))*fact;
		sa += fact;
	}

	X(0) = cx/(3.*sa);
	X(1) = cy/(3.*sa);
	X(2) = ver[0](2);

	// cout << "X: " << X.transpose() << endl;

	X = Qf._transformVector(X);
	// cout << "X: " << X.transpose() << endl;
	P0.resize(nv);
	for (size_t i=0; i<nv; ++i)
	{
		P0[i] = Qfi._transformVector(ver0[i]-X);
	}
}

inline void DEM_PARTICLE::SetSphere(double r)
{
    ShapeType= 1;
	R	= r;
	Vol = 4./3.*M_PI*R*R*R;
	M	= Rho*Vol;
	I(0)	= 0.4*M*R*R;
	I(1)	= I(0);
	I(2)	= I(0);

	P0.resize(1);
	P0[0] << 0., 0., 0.;

	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	BoxL << R+1, R+1, R+1;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));

	MetaP0.push_back(P0[0]);
	MetaK.resize(1);
	MetaK(0) = R*R;
}
// only used for 2d
inline void DEM_PARTICLE::SetDisk2D(double r)
{
    ShapeType= 1;
	R	= r;
	Vol = M_PI*R*R;
	M	= Rho*Vol;
	I(0)	= 0.25*M*R*R;
	I(1)	= I(0);
	I(2)	= 2.*I(0);

	P0.resize(1);
	P0[0] << 0., 0., 0.;

	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	BoxL << R+1, R+1, R+1;

	Nfe = 0;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));

	MetaP0.push_back(P0[0]);
	MetaK.resize(1);
	MetaK(0) = R*R;
}

inline void DEM_PARTICLE::SetCuboid(double lx, double ly, double lz)
{
    ShapeType= 2;
	R	= 0.501*sqrt(lx*lx+ly*ly+lz*lz);
	Vol = lx*ly*lz;
	M	= Rho*Vol;
	I(0)	= M*(ly*ly+lz*lz)/12.;
	I(1)	= M*(lx*lx+lz*lz)/12.;
	I(2)	= M*(lx*lx+ly*ly)/12.;

	P0.resize(8);
	P0[0] << -0.5*lx, -0.5*ly, -0.5*lz;
	P0[1] <<  0.5*lx, -0.5*ly, -0.5*lz;
	P0[2] <<  0.5*lx,  0.5*ly, -0.5*lz;
	P0[3] << -0.5*lx,  0.5*ly, -0.5*lz;
	P0[4] << -0.5*lx, -0.5*ly,  0.5*lz;
	P0[5] <<  0.5*lx, -0.5*ly,  0.5*lz;
	P0[6] <<  0.5*lx,  0.5*ly,  0.5*lz;
	P0[7] << -0.5*lx,  0.5*ly,  0.5*lz;

	VectorXi face(4);
	
	face << 0, 3, 2, 1;
	Faces.push_back(face);
	face << 4, 5, 6, 7;
	Faces.push_back(face);
	face << 1, 5, 4, 0;
	Faces.push_back(face);
	face << 2, 3, 7, 6;
	Faces.push_back(face);
	face << 1, 2, 6, 5;
	Faces.push_back(face);
	face << 0, 4, 7, 3;
	Faces.push_back(face);

	Nfe = 30;

	Q0.w() = 1;
	Q0.vec() << 0.,0.,0.;

	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	BoxL << R+1, R+1, R+1;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetTetrahedron(vector<Vector3d> ver)
{
    ShapeType= 2;
	Vol	= abs((ver[0]-ver[3]).dot((ver[1]-ver[3]).cross(ver[2]-ver[3])))/6.;
	M	= Rho*Vol;
	X 	= (ver[0]+ver[1]+ver[2]+ver[3])/4.;
	// add Faces
	VectorXi face(3);
	face << 0, 1, 2;
	Faces.push_back(face);
	face << 1, 2, 3;
	Faces.push_back(face);
	face << 2, 3, 0;
	Faces.push_back(face);
	face << 0, 3, 1;
	Faces.push_back(face);
	RerangeFaceElementOrder();

	Matrix3d Ip;	// Inertia tensor
	CalTriangleMeshProperties(Faces, ver, Rho, Vol, M, X, Ip);
	// use eigen vectors to find the rotation from lab frame to object frame
	SelfAdjointEigenSolver<Matrix3d> sol (Ip);
	Matrix3d em = sol.eigenvectors();
	if ((em.col(0).cross(em.col(1))).dot(em.col(2))<0)	em.col(2) *= -1;
	Quaterniond q0(em);

	Q0 = q0;
	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();
	// move P0 to make sure mass centre is the origin
	P0.resize(ver.size());
	double maxr = 0.;
	for (size_t i=0; i<P0.size(); ++i)
	{
		P0[i] = ver[i]-X;
		P0[i] = q0.inverse()._transformVector(P0[i]);
		double r = P0[i].norm();
		if (r>maxr)	maxr = r;
	}
	R = 1.1*maxr;

	Nfe = 16;

	BoxL << R+1, R+1, R+1;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetMetaball(double dis, vector<Vector3d>& metaP, VectorXd& metaK)
{
	ShapeType = 4;
	MetaK = metaK;

	vector<Vector3d> ver(0);
	CalSurfaceMesh(metaP, MetaK, dis, ver, Faces);
	Matrix3d Ip;	// Inertia tensor
	CalTriangleMeshProperties(Faces, ver, Rho, Vol, M, X, Ip);
	// use eigen vectors to find the rotation from lab frame to object frame
	SelfAdjointEigenSolver<Matrix3d> sol (Ip);
	Matrix3d em = sol.eigenvectors();
	I = sol.eigenvalues();
	if ((em.col(0).cross(em.col(1))).dot(em.col(2))<0)	em.col(2) *= -1;

	// connectivity for points
	PNei.resize(ver.size());
	Parea.resize(ver.size(),0.);

	for (size_t i=0; i<Faces.size(); ++i)
	{
		int p0 = Faces[i](0);
		int p1 = Faces[i](1);
		int p2 = Faces[i](2);
		// area for each point
		double areai = 1./6.*((ver[p1]-ver[p0]).cross(ver[p2]-ver[p0])).norm();
		Parea[p0] += areai;
		Parea[p1] += areai;
		Parea[p2] += areai;

		if (p0>p1)
		{
			PNei[p0].push_back(p1);
			PNei[p1].push_back(p0);
		}
		if (p1>p2)
		{
			PNei[p1].push_back(p2);
			PNei[p2].push_back(p1);
		}
		if (p2>p0)
		{
			PNei[p2].push_back(p0);
			PNei[p0].push_back(p2);
		}
	}

	Quaterniond q0(em);
	// move P0 to make sure mass centre is the origin and rotation P0 to object frame
	P0.resize(ver.size());
	double maxr = 0.;
	for (size_t i=0; i<ver.size(); ++i)
	{
		// correct errors between surface point and metaball function
		double c = CalMetaC(metaP, metaK, ver[i]);
		if (abs(c-1.)>1.e-12)
		{
			double ci;
			Vector3d fi;
			CalMetaCF(metaP, metaK, ver[i], ci, fi);
			double q = (1.-ci)/fi.squaredNorm();
			ver[i] += q*fi;
		}
		P0[i] = ver[i]-X;
		P0[i] = q0.inverse()._transformVector(P0[i]);
		double r = P0[i].norm();
		if (r>maxr)	maxr = r;
	}
	R = 1.1*maxr;
	MetaP0.resize(metaP.size());
	for (size_t i=0; i<metaP.size(); ++i)
	{
		MetaP0[i] = metaP[i]-X;
		MetaP0[i] = q0.inverse()._transformVector(MetaP0[i]);
	}

	Nfe = 4*Faces.size();

	BoxL << R+1, R+1, R+1;

	Q0 = q0;
	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

// warning: polygon2D cannot have internal holes
inline void DEM_PARTICLE::SetPolygon2D(vector<Vector3d> ver)
{
	ShapeType= 7;
	size_t nv = ver.size();	// number of vertices
	double cx = 0.; double cy = 0.;
	double jx = 0.; double jy = 0.;
	double sa = 0.;
	for (size_t i=0; i<nv; ++i)
	{
		size_t j = (i+1)%nv;
		double fact = (ver[i](0)*ver[j](1) - ver[j](0)*ver[i](1));
		cx += (ver[i](0)+ver[j](0))*fact;
		cy += (ver[i](1)+ver[j](1))*fact;
		jx += (ver[i](1)*ver[i](1) + ver[i](1)*ver[j](1) + ver[j](1)*ver[j](1))*fact;
		jy += (ver[i](0)*ver[i](0) + ver[i](0)*ver[j](0) + ver[j](0)*ver[j](0))*fact;
		sa += fact;
	}

	X(0) = cx/(3.*sa);
	X(1) = cy/(3.*sa);
	X(2) = 0.;

	Vol = 0.5*sa;
	M	= Rho*Vol;

	I(0)	= Rho*jx/12.;	// not correct but is ok for 2d
	I(1)	= Rho*jy/12.;	// not correct but is ok for 2d
	I(2)	= I(0)+I(1) - M*X.squaredNorm();

	Vector3d n0 (0.,0.,1.);
	Quaterniond q0;
	q0.setFromTwoVectors(n0,n0);	
	Q0 = q0;

	P0.resize(nv);
	VectorXi face(nv);
	double maxDis = 0.;
	for (size_t i=0; i<nv; ++i)
	{
		P0[i] = ver[i]-X;
		if (P0[i].norm()>maxDis)	maxDis = P0[i].norm();
		size_t j = (i+1)%nv;
		Vector2i edge (i, j);
		Edges.push_back(edge);
		face(i) = i;
	}
	Faces.push_back(face);
	Nfe = nv+1;

	R = maxDis;

	BoxL << R, R, 0.;

	Q0.w() = 1;
	Q0.vec() << 0.,0.,0.;

	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetCylinder(double h, double r, Vector3d n)
{
    ShapeType= 3;
	R	= sqrt(0.25*h*h + r*r);
	Vol = h*M_PI*r*r;
	M	= Rho*Vol;
	I(0)	= M*(3.*r*r + h*h)/12.;
	I(1)	= I(0);
	I(2)	= 0.5*M*r*r;

	Nfe = 0;

	size_t num = 50;
	P0.resize(2*(num+1));
	// P.resize(2*(num+1));
	double ang = 2.*M_PI/num;
	for (size_t i=0; i<num; ++i)
	{
		P0[2*i] << cos(i*ang)*r, sin(i*ang)*r, 0.5*h;
		P0[2*i+1] << cos(i*ang)*r, sin(i*ang)*r, -0.5*h;

		VectorXi face(4);
		face << 2*i, 2*((i+1)%num), 2*((i+1)%num)+1, 2*i+1;
		Faces.push_back(face);
		Nfe += 5;
	}
	size_t indT = P0.size()-2;
	size_t indB = P0.size()-1;

	P0[indT] << 0., 0., 0.5*h;
	P0[indB] << 0., 0., -0.5*h;

	Vector3d nc = n.normalized();
	Vector3d n0 (0.,0.,1.);
	Quaterniond q0;
	q0.setFromTwoVectors(n0,nc);
	Q0 = q0;
	Qf = Q0*Q;
	Qfi = Qf.inverse();

	for (size_t i=0; i<num; ++i)
	{
		VectorXi face(3);
		face << indT, 2*i, 2*((i+1)%num);
		Faces.push_back(face);
		Nfe += 4;
	}
	for (size_t i=0; i<num; ++i)
	{
		VectorXi face(3);
		face << indB, 2*((i+1)%num)+1, 2*i+1;
		Faces.push_back(face);
		Nfe += 4;
	}
	RerangeFaceElementOrder();

	BoxL << R+1, R+1, R+1;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetTriangle2D(vector<Vector3d> ver)
{
	ShapeType= 6;
	Vol = ((ver[1]-ver[0]).cross(ver[2]-ver[0])).norm()/2.;
	M	= Rho*Vol;
	X 	= (ver[0]+ver[1]+ver[2])/3.;
	double a2 = (ver[1]-ver[0]).squaredNorm();
	double b2 = (ver[2]-ver[0]).squaredNorm();
	double c2 = (ver[1]-ver[2]).squaredNorm();
	I(2)	= M*(a2+b2+c2)/36.;
	I(1)	= I(2);		// not correct but is ok for 2d
	I(0)	= I(2);		// not correct but is ok for 2d

	// P = ver;
	P0.resize(ver.size());
	for (size_t i=0; i<ver.size(); ++i)		P0[i] = ver[i]-X;
	// add Edges
	Vector2i edge;
	edge << 0, 1;
	Edges.push_back(edge);
	edge << 1, 2;
	Edges.push_back(edge);
	edge << 2, 0;
	Edges.push_back(edge);
	// add Faces
	VectorXi face(3);
	face << 0, 1, 2;
	Faces.push_back(face);
	Nfe = 4;

	R = max(max((ver[0]-X).norm(), (ver[1]-X).norm()), (ver[2]-X).norm());

	Q0.w() = 1;
	Q0.vec() << 0.,0.,0.;

	Qf = Q0*Q;	// final rotation (from object frame to lab frame)
	Qfi = Qf.inverse();

	BoxL << R+1, R+1, 0.;

	Max(0) 	= (int) (X(0)+BoxL(0));
	Max(1) 	= (int) (X(1)+BoxL(1));
	Max(2) 	= (int) (X(2)+BoxL(2));

	Min(0) 	= (int) (X(0)-BoxL(0));
	Min(1) 	= (int) (X(1)-BoxL(1));
	Min(2) 	= (int) (X(2)-BoxL(2));
}

inline void DEM_PARTICLE::SetFromINP(string fname)
{
	cout << "Start reading from Abaqus inp: " << fname << endl;
	ifstream in(fname);
	if (!in)
	{
		std::cout << "Could not open file: " << fname << std::endl;
		abort();
	}

	string dimTag = "ELSET=Surface1";
	string dimTagSurf = "ELSET=Line";
	// if (D==3)
	{
		dimTag = "ELSET=Volume1";
		dimTagSurf = "ELSET=Surface";
	}

	int numSurfNode = 0;

	string line;
	double i, x, y, z;
	int type = -1;

	vector<size_t>	pid;
	vector<Vector3d> pos;
	vector<VectorXi> ele;
	vector<VectorXi> surf;

	int eleType = -1;

	while ( getline( in, line ) )                          
	{
		if (line[0]!='*')
		{
			if (type==0)	// read node
			{
		        vector<double> subStr;
		        stringstream sStream(line);
				while(sStream.good())
				{
					string subStri;
					getline(sStream, subStri, ','); //get first string delimited by comma
					double numi = stod(subStri);
					subStr.push_back(numi);
				}

				i = ((size_t) subStr[0])-1;
				x = subStr[1];
				y = subStr[2];
				z = subStr[3];

				pid.push_back(i);
				pos.push_back({x,y,z});
			}	
			else if (type==1)	// read face
			{
		        vector<int> subStr;
		        stringstream sStream(line);
				while(sStream.good())
				{
					string subStri;
					getline(sStream, subStri, ','); //get first string delimited by comma
					int numi = stoi(subStri);
					subStr.push_back(numi);
				}
				int numNPE = subStr.size()-1;	// number of node per element
				VectorXi elem(numNPE);
				if (eleType==-1)	eleType = numNPE;
				else if (numNPE!=eleType)
				{
					cout << "number of node per element: " << numNPE << endl;
					cout << "eleType: " << eleType << endl;
					cout << "don't support mixed element type yet." << endl;
					abort();
				}
				for (size_t m=1; m<subStr.size(); ++m)
				{
					elem[m-1] = subStr[m]-1;
				}
				ele.push_back(elem);
			}
			else if (type==2)
			{
		        vector<int> subStr;
		        stringstream sStream(line);
				while(sStream.good())
				{
					string subStri;
					getline(sStream, subStri, ','); //get first string delimited by comma
					int numi = stoi(subStri);
					subStr.push_back(numi);
				}
				VectorXi surfm(numSurfNode);
				for (size_t m=1; m<subStr.size(); ++m)
				{
					surfm[m-1] = subStr[m]-1;
				}
				surf.push_back(surfm);
			}
		}
		else	// for infor line
		{
			type = -1;		// reset flag
		  	if (line[1]=='N' && line[2]=='O')	// for node
		  	{
		  		type = 0;
		  	}
		  	else if (line[1]!='*' && line[1]!='H')
		  	{
		        vector<string> subStr;
		        stringstream sStream(line);
				while(sStream.good()) 
				{
					string subStri;
					getline(sStream, subStri, ','); //get first string delimited by comma
					subStr.push_back(subStri);
				}
				if (subStr[0]=="*ELEMENT")
				{
					string tag0 = subStr[2];
					tag0.erase(0, 1);
					tag0.resize(dimTag.size());

					if (tag0 == dimTag)	type = 1;		// only one surface/volume
					else
					{
						string tag1 = subStr[2];
						tag1.erase(0, 1);
						tag1.resize(dimTagSurf.size());
						if (tag1 == dimTagSurf)
						{
							type = 2;
							string tag2 = subStr[1];
							numSurfNode = (int) tag2.back()-48;
						}
					}
				}
		  	}
		}
	}
	in.close();

	unordered_map<size_t, size_t> pidMap(0);
	for (size_t i=0; i<pid.size(); ++i)
	{
		pidMap.insert({pid[i],i});
	}

	for (size_t i=0; i<surf.size(); ++i)
	{
		for (int j=0; j<surf[i].size(); ++j)
		{
			size_t id = pidMap[surf[i](j)];
			surf[i](j) = id;
		}
	}


	X.setZero();

	Faces.clear();
	P0.clear();

	Faces = surf;
	P0 = pos;

	Nfe = 4*Faces.size();


	cout << "pos.size(): " << pos.size() << endl;
	cout << "surf.size(): " << surf.size() << endl;

	cout << "Reading mesh finshed." << endl;
}

#endif