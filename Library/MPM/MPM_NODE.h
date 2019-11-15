class MPM_NODE
{
public:
	MPM_NODE();
	MPM_NODE(size_t level, const Vector3d& x);
	void Reset();
	void ResetwithDEM(size_t nproc);
	void CombineDPs();
	void NonSlippingBC();
	void SlippingBC(Vector3d& norm);
	void FrictionBC(double dt, Vector3d& n);

	bool 							Actived;					// Flag of actived node

	size_t 							Type;						// 0 for nodes, 1 for gauss quadrature points
	size_t 							ID; 				    	// Index of gird in the list 
	size_t 							Level; 				    	// Level of the node, 0 for basic node, 1 for level 1 refined node

	double							M;							// Mass
	double							Vol;						// Volume
	double 							Mu;							// Friction coefficient

	Vector3d						X;							// Position
	Vector3d						V;							// Velocity
	Vector3d						Mv;							// Momentum
	Vector3d						F;							// Total force
	Vector3d						Fi;							// Internal force
	Vector3d						Fe;							// External force

	Matrix3d						S;							// Stress

	vector<size_t>					GPs;						// Gauss quadrature points
	vector<double> 					Ng;							// Shape function
	vector<Vector3d>				GNg;						// Gradient of shape function
	vector<size_t>					MPs;						// Material points
	vector<size_t>					DPs;						// DEM particles
	vector<vector<size_t>>			DPs_proc;					// DEM particles for openMP to avoid race condition
	vector<size_t> 					BCTypes;					// Boundary condition type
	vector<Vector3d> 				Norms;						// Normal direction for boundaries

	// unordered_map<size_t, bool> 	PMap;						// Map for all particles which have influence on this node
};

inline MPM_NODE::MPM_NODE()
{
	Actived = false;
	Type 	= 0;
	ID 		= 0;
	Level 	= 0;
	M 		= 0.;
	Vol 	= 1.;
	X 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	GPs.resize(0);
	MPs.resize(0);
	DPs.resize(0);
	DPs_proc.resize(0);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline MPM_NODE::MPM_NODE(size_t level, const Vector3d& x)
{
	Actived = false;
	Type 	= 0;
	ID 		= 0;
	Level 	= level;
	M 		= 0.;
	Vol 	= 1.;
	X 		= x;
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
	Fi 		= Vector3d::Zero();
	Fe 		= Vector3d::Zero();
	GPs.resize(0);
	MPs.resize(0);
	DPs.resize(0);
	DPs_proc.resize(0);
	BCTypes.resize(0);
	Norms.resize(0);
}

inline void MPM_NODE::Reset()
{
	Actived = false;
	M = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
	Fi.setZero();
	Fe.setZero();
	S.setZero();
	// MPs.resize(0);
	MPs.clear();
}

inline void MPM_NODE::ResetwithDEM(size_t nproc)
{
	Actived = false;
	M = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
	Fi.setZero();
	Fe.setZero();
	S.setZero();
	MPs.clear();
	DPs.clear();
	// DPs_proc.resize(nproc);
	for (size_t i=0; i<DPs_proc.size(); ++i)
	{
		DPs_proc[i].clear();
	}
}

inline void MPM_NODE::CombineDPs()
{
	// no need to sort the list since DEM is parallized base on particles
	for (size_t i=0; i<DPs_proc.size(); ++i)
	{
		DPs.insert( DPs.end(), DPs_proc[i].begin(), DPs_proc[i].end() );
	}
	// sort( DPs.begin(), DPs.end() );
	// DPs.erase(unique(DPs.begin(), DPs.end()), DPs.end());
}

inline void MPM_NODE::NonSlippingBC()
{
	F.setZero();
	Mv.setZero();
}

inline void MPM_NODE::SlippingBC(Vector3d& n)
{	
	// make sure it's compressing the bc then set norm force and momentum to zero
	if (Mv.dot(n)>0.)
	{
		F  = F-F.dot(n)*n;
		Mv = Mv-Mv.dot(n)*n;
	}
}

inline void MPM_NODE::FrictionBC(double dt, Vector3d& n)
{
	double fnNorm = F.dot(n);							// Normal force magnitude
	if (fnNorm>0.)
	{
		Vector3d vt0 = V-V.dot(n)*n;					// Tangential velocity without friction (5.9)
		Vector3d mvn0 = Mv.dot(n)*n;					// Normal momentum
		Vector3d mvt0 = Mv-mvn0;						// Tangential momentum
		Vector3d t = mvt0.normalized();					// Tangential vector
		Vector3d ff = -Sign(vt0.dot(t))*Mu*fnNorm*t;	// Friction force (5.10)
		Vector3d mvt1 = mvt0+ff*dt;						// (5.10)
		if (mvt0.dot(mvt1)>0.)
		{
			F = F.dot(t)*t+ff;
			Mv = mvt1;
		}
		else
		{
			F.setZero();
			Mv.setZero();
		}
	}
}