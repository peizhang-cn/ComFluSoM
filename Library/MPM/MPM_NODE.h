class MPM_NODE
{
public:
	MPM_NODE();
	MPM_NODE(size_t level, const Vector3d& x);
	void Reset();
	void NonSlippingBC();
	void SlippingBC(Vector3d& norm);
	void FrictionBC(double dt, Vector3d& n);

	size_t 							ID; 				    	// Index of gird in the list 
	size_t 							Level; 				    	// Level of the node, 0 for basic node, 1 for level 1 refined node

	double							M;							// Mass
	double 							Mu;							// Friction coefficient

	Vector3d						X;							// Position
	Vector3d						V;							// Velocity
	Vector3d						Mv;							// Momentum
	Vector3d						F;							// Total force

	vector<size_t> 					BCTypes;					// Boundary condition type
	vector<Vector3d> 				Norms;						// Normal direction for boundaries

	unordered_map<size_t, bool> 	PMap;						// Map for all particles which have influence on this node
};

inline MPM_NODE::MPM_NODE()
{
	ID 		= 0;
	Level 	= 0;
	M 		= 0.;
	X 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();

	BCTypes.resize(0);
	Norms.resize(0);
}

inline MPM_NODE::MPM_NODE(size_t level, const Vector3d& x)
{
	ID 		= 0;
	Level 	= level;
	M 		= 0.;
	X 		= x;
	V 		= Vector3d::Zero();
	Mv 		= Vector3d::Zero();
	F 		= Vector3d::Zero();
}

inline void MPM_NODE::Reset()
{
	M = 0.;
	V.setZero();
	Mv.setZero();
	F.setZero();
}

inline void MPM_NODE::NonSlippingBC()
{
	F.setZero();
	Mv.setZero();
}

inline void MPM_NODE::SlippingBC(Vector3d& n)
{	
	F  = F-F.dot(n)*n;
	Mv = Mv-Mv.dot(n)*n; 
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