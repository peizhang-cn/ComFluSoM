// GJK algorithm for collision detection and shortest distance calculatation between convex polyhedron

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

// Support function
void Support(vector<Vector3d>& veri, vector<Vector3d>& verj, Vector3d& n, size_t i, size_t j, Vector3d& p)
{
	i = FarthestPoint(veri, n);
	j = FarthestPoint(verj, -n);
	p = veri[i]-verj[j];
}

void FindClosestPoints3D(vector<Vector3d>& veri, vector<Vector3d>& verj, Vector3d& n)
{
	size_t i0, j0, i1, j1, i2, j2, i3, j3;
	Vector3d p0, p1, p2, p3;
	Support(veri, verj, n, i0, j0, p0);
	n = -n;
	Support(veri, verj, n, i1, j1, p1);
	n = -p1;
	Support(veri, verj, n, i2, j2, p2);
	if (p2==p0)
	{
		cout << "p2=p0, cannot form trangle." << endl;
		// abort();
	}
	cout << p0.transpose() << endl;
	cout << p1.transpose() << endl;
	cout << p2.transpose() << endl;

	n = (p1-p0).cross(p2-p0);
	if (-p0.dot(n)<0.)	n *= -1;
	n.normalize();
	cout << n.transpose() << endl;
	Support(veri, verj, n, i3, j3, p3);
	cout << "p3: " << p3.transpose() << endl;
}