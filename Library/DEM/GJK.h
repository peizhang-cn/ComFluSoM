// GJK algorithm for shortest distance calculatation between non-overlapping convex polyhedron
// NOTE: This is not a general GJK since only NON-OVERLAPPING cases are considered here (07 Feb 2020).

struct CSOPoint
{
	Vector3d X;
	size_t I;
	size_t J;
};

// Return the index of farthest point in n direction
size_t FarthestPoint(/*vertices*/vector<Vector3d> ver, /*dirction*/Vector3d n)
{
	size_t ind = 0;
	double maxDot = ver[0].dot(n);
	for (size_t i=1; i<ver.size(); ++i)
	{
		double dot = ver[i].dot(n);
		if (dot>maxDot)
		{
			ind = i;
			maxDot = dot;
		}
	}
	return ind;
}

CSOPoint Support(vector<Vector3d>& veri, vector<Vector3d>& verj, Vector3d& n)
{
	CSOPoint p;
	p.I = FarthestPoint(veri, n);
	p.J = FarthestPoint(verj, -n);
	p.X = veri[p.I]-verj[p.J];
	return p;
}

// Find the closest point from a line segment to origin
Vector3d LineClosestPointToOrigin(Vector3d& A, Vector3d& B)
{
	Vector3d AB = B-A;
	double k = -A.dot(AB)/AB.squaredNorm();
	Vector3d p = A;
	if (k>1.)		p = B;
	else if (k>0.)	p += k*AB;
	return p;
}

Vector3d TriangleClosestPointToOrigin(/*the newest point*/Vector3d& A, Vector3d& B, Vector3d& C)
{
	Vector3d AC = C-A;
	Vector3d AB = B-A;
	Vector3d ABC = AB.cross(AC);
	Vector3d p;

	if (A.dot(ABC.cross(AC))<0.)		p = LineClosestPointToOrigin(A, C);
	else if (A.dot(AB.cross(ABC))<0.)	p = LineClosestPointToOrigin(A, B);
	else								p = ABC.dot(A)/ABC.squaredNorm()*ABC;
	return p;	
}

// void TetrahedronDirectionToOrigin(/*the newest point*/Vector3d& A, Vector3d& B, Vector3d& C, Vector3d& D, Vector3d& n, size_t& ind0)
// {
// 	Vector3d BCD = (B-C).cross(D-C);
// 	cout << "dot: " << ((B-C).cross(D-C)).dot(A-C) << endl;
// 	if (((B-C).cross(D-C)).dot(A-C) == 0.)
// 	{
// 		cout << "intointointointointointointointointointointointointointointointointointointo" << endl;
// 		Vector3d c0 = TriangleClosestPointToOrigin(A, B, C);
// 		Vector3d c1 = TriangleClosestPointToOrigin(A, B, D);
// 		Vector3d c2 = TriangleClosestPointToOrigin(A, C, D);
// 		double cn0 = c0.norm(); double cn1 = c1.norm(); double cn2 = c2.norm();
// 		double mind = min(cn0, min(cn1,cn2));
// 		if (cn0==mind)
// 		{
// 			n = -c0;
// 			ind0 = 0;
// 		}
// 		else if (cn1==mind)
// 		{
// 			n = -c1;
// 			ind0 = 1;
// 		}
// 		else if (cn2==mind)
// 		{
// 			n = -c2;
// 			ind0 = 2;
// 		}
// 	}
// 	else
// 	{
// 		// find norm direction (towards to outside) of face ABC
// 		Vector3d ABC = (B-A).cross(C-A);
// 		if ((D-A).dot(ABC)>0.)	ABC *= -1;
// 		if (A.dot(ABC)<0.)
// 		{
// 			n = -TriangleClosestPointToOrigin(A, B, C);
// 			ind0 = 0;
// 		}
// 		else
// 		{
// 			Vector3d ACD = (C-A).cross(D-A);
// 			if ((B-A).dot(ACD)>0.)	ACD *= -1;
// 			if (A.dot(ACD)<0.)
// 			{
// 				n = -TriangleClosestPointToOrigin(A, C, D);
// 				ind0 = 2;
// 			}
// 			else
// 			{
// 				n = -TriangleClosestPointToOrigin(A, B, D);
// 				ind0 = 1;
// 			}
// 		}
// 	}
// }

void TetrahedronDirectionToOrigin(/*the newest point*/Vector3d& A, Vector3d& B, Vector3d& C, Vector3d& D, Vector3d& n, size_t& ind0)
{
		// cout << "intointointointointointointointointointointointointointointointointointointo" << endl;
		Vector3d c0 = TriangleClosestPointToOrigin(A, B, C);
		Vector3d c1 = TriangleClosestPointToOrigin(A, B, D);
		Vector3d c2 = TriangleClosestPointToOrigin(A, C, D);
		double cn0 = c0.norm(); double cn1 = c1.norm(); double cn2 = c2.norm();
		double mind = min(cn0, min(cn1,cn2));
		if (cn0==mind)
		{
			n = -c0;
			ind0 = 0;
		}
		else if (cn1==mind)
		{
			n = -c1;
			ind0 = 1;
		}
		else if (cn2==mind)
		{
			n = -c2;
			ind0 = 2;
		}
}

bool CheckSimplex(vector<CSOPoint> s)
{
	bool tst = false;
	for (size_t i=0; i<s.size()-1; ++i)
	{
		if (s[s.size()-1].X==s[i].X)
		{
			tst = true;
			break;
		}
	}
	return tst;
}

void showFunc(vector<Vector3d>& veri, vector<Vector3d>& verj, /*initial direction*/Vector3d& n, Vector3d& pi, Vector3d& pj)
{
	vector<CSOPoint> simplex;
	// Add first point to simplex
	simplex.push_back(Support(veri, verj, n));
	// Add second point to simplex
	n = -simplex[0].X;
	simplex.push_back(Support(veri, verj, n));
	// Add third point to simplex
	n = -LineClosestPointToOrigin(simplex[0].X, simplex[1].X);
	cout << "n line: " << n.transpose() << endl;
	simplex.push_back(Support(veri, verj, n));

	if (CheckSimplex(simplex))
	{
		double lamda = (-n-simplex[1].X).norm()/(simplex[0].X-simplex[1].X).norm();
		cout << "lamda: " << lamda << endl;
		cout << "simplex[0].X: " << simplex[0].X.transpose() << endl;
		cout << "simplex[1].X: " << simplex[1].X.transpose() << endl;
		cout << "-n: " << -n.transpose() << endl;

		cout << "cal n: " << (lamda*simplex[0].X+(1.-lamda)*simplex[1].X).transpose() << endl;

		cout << veri[simplex[0].I].transpose() << endl;
		cout << veri[simplex[1].I].transpose() << endl;
		cout << "--------------------" << endl;
		cout << verj[simplex[0].J].transpose() << endl;
		cout << verj[simplex[1].J].transpose() << endl;

		pi = lamda*veri[simplex[0].I] + (1.-lamda)*veri[simplex[1].I];
		pj = lamda*verj[simplex[0].J] + (1.-lamda)*verj[simplex[1].J];
		// simplex.pop_back();
		// Vector3d ni = -n;
		// simplex.push_back(Support(veri, verj, ni));
		// if (CheckSimplex(simplex))
		// {
		// 	cout << "still cannot form Triangle" << endl;
		// }
	}
	else
	{
		Vector3d A = simplex[2].X;
		Vector3d B = simplex[1].X;
		Vector3d C = simplex[0].X;

		cout << "A: " << A << endl;
		cout << "B: " << B << endl;
		cout << "C: " << C << endl;

		n = -TriangleClosestPointToOrigin(A, B, C);
		cout << "n: " << n.transpose() << endl;
		simplex.push_back(Support(veri, verj, n));

		cout << "D: " << simplex[3].X.transpose() << endl;

		size_t ind0;
		size_t count = 0;
		while (!CheckSimplex(simplex))
		{
			cout << "befrore simplex[3].X: " << simplex[3].X.transpose() << endl;
			cout << "befrore simplex[2].X: " << simplex[2].X.transpose() << endl;
			cout << "befrore simplex[1].X: " << simplex[1].X.transpose() << endl;
			cout << "befrore simplex[0].X: " << simplex[0].X.transpose() << endl;
			cout << "----------------------------------------------" << endl;
			Vector3d nb = n;
			TetrahedronDirectionToOrigin(simplex[3].X, simplex[2].X, simplex[1].X, simplex[0].X, n, ind0);
			cout << "ind0= " << ind0 << endl;
			simplex.erase(simplex.begin()+ind0);
			simplex.push_back(Support(veri, verj, n));
			cout << "simplex[3].X: " << simplex[3].X.transpose() << endl;
			cout << "simplex[2].X: " << simplex[2].X.transpose() << endl;
			cout << "simplex[1].X: " << simplex[1].X.transpose() << endl;
			cout << "simplex[0].X: " << simplex[0].X.transpose() << endl;

			cout << "n.norm(): " << n.norm() << endl;

			cout << "ang: " << acos(nb.dot(n)/nb.norm()/n.norm())*180/M_PI << endl;
			cout << "============================================" << endl;
			count++;
			if (count>100)
			{
				cout << "loop 100 times and still not find GJK solution!" << endl;
				cout << "veri[0]: "<< veri[0].transpose() << endl;
				cout << "verj[0]: "<< verj[0].transpose() << endl;
				abort();
			}	
		}

		// for (size_t i=0; i<simplex.size(); ++i)
		// {
		// 	cout << simplex[i].X.transpose() << endl;
		// }

		Vector3d cp = TriangleClosestPointToOrigin(simplex[2].X, simplex[1].X, simplex[0].X);

		Matrix3d mat;
		mat.col(0) = simplex[0].X;
		mat.col(1) = simplex[1].X;
		mat.col(2) = simplex[2].X;
		Vector3d coef = mat.inverse()*cp;
		cout << "coef: " << coef.transpose() << endl;

		pi = coef(0)*veri[simplex[0].I] + coef(1)*veri[simplex[1].I] + coef(2)*veri[simplex[2].I];
		pj = coef(0)*verj[simplex[0].J] + coef(1)*verj[simplex[1].J] + coef(2)*verj[simplex[2].J];

		cout << veri[simplex[0].I].transpose() << endl;
		cout << veri[simplex[1].I].transpose() << endl;
		cout << veri[simplex[2].I].transpose() << endl;
		cout << "--------------------" << endl;
		cout << verj[simplex[0].J].transpose() << endl;
		cout << verj[simplex[1].J].transpose() << endl;
		cout << verj[simplex[2].J].transpose() << endl;
	}

	cout << "cpi: " << pi.transpose() << endl;
	cout << "cpj: " << pj.transpose() << endl;
}

void FindClosestPoints3D(vector<Vector3d>& veri, vector<Vector3d>& verj, /*initial direction*/Vector3d& n, Vector3d& pi, Vector3d& pj)
{
	Vector3d n0 = n;
	vector<CSOPoint> simplex;
	// Add first point to simplex
	simplex.push_back(Support(veri, verj, n));
	// Add second point to simplex
	n = -simplex[0].X;
	simplex.push_back(Support(veri, verj, n));
	// Add third point to simplex
	n = -LineClosestPointToOrigin(simplex[0].X, simplex[1].X);;
	simplex.push_back(Support(veri, verj, n));

	Vector3d A = simplex[2].X;
	Vector3d B = simplex[1].X;
	Vector3d C = simplex[0].X;

	n = -TriangleClosestPointToOrigin(A, B, C);
	simplex.push_back(Support(veri, verj, n));

	size_t ind0;
	size_t count = 0;
	while (!CheckSimplex(simplex))
	{
		TetrahedronDirectionToOrigin(simplex[3].X, simplex[2].X, simplex[1].X, simplex[0].X, n, ind0);
		simplex.erase(simplex.begin()+ind0);
		double disb = n.norm();
		simplex.push_back(Support(veri, verj, n));
		double dis = n.norm();
		count++;
		if (count>100)
		{
			cout << "loop 100 times and still not find GJK solution!" << endl;
			cout << "veri[0]: "<< veri[0].transpose() << endl;
			cout << "verj[0]: "<< verj[0].transpose() << endl;
			Vector3d pi0, pj0;
			showFunc(veri, verj, n0, pi0, pj0);

			abort();
		}
		if (dis<disb)
		{
			cout << "dis<disb" << endl;
		}	
	}

	// for (size_t i=0; i<simplex.size(); ++i)
	// {
	// 	cout << simplex[i].X.transpose() << endl;
	// }

	Vector3d cp = TriangleClosestPointToOrigin(simplex[2].X, simplex[1].X, simplex[0].X);

	Matrix3d mat;
	mat.col(0) = simplex[0].X;
	mat.col(1) = simplex[1].X;
	mat.col(2) = simplex[2].X;
	Vector3d coef = mat.inverse()*cp;
	cout << "coef: " << coef.transpose() << endl;

	pi = coef(0)*veri[simplex[0].I] + coef(1)*veri[simplex[1].I] + coef(2)*veri[simplex[2].I];
	pj = coef(0)*verj[simplex[0].J] + coef(1)*verj[simplex[1].J] + coef(2)*verj[simplex[2].J];


	cout << veri[simplex[0].I].transpose() << endl;
	cout << veri[simplex[1].I].transpose() << endl;
	cout << veri[simplex[2].I].transpose() << endl;
	cout << "--------------------" << endl;
	cout << verj[simplex[0].J].transpose() << endl;
	cout << verj[simplex[1].J].transpose() << endl;
	cout << verj[simplex[2].J].transpose() << endl;

	cout << "cpi: " << pi.transpose() << endl;
	cout << "cpj: " << pj.transpose() << endl;
}