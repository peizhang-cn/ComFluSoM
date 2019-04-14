class RWM_PARTICLE
{
public:
	RWM_PARTICLE();
	RWM_PARTICLE(int type, const Vector3d& x, double m);

    int         				Type;                       // Type of particle, for -1 is freely moved particles or -2 is boundary particles.
	int         				ID; 				    	// Index of particle in the list 
	int         				Tag;				    	// Tag of particle

	double    					M;				            // Mass

	Vector3d    				X;				            // Position
	Vector3d    				Xb;				            // Position one time step before

    bool        				Removed;                	// Flag for removed particles
};

inline RWM_PARTICLE::RWM_PARTICLE()
{
    Type	= -1;
	ID		= 0;
	Tag		= 0;
	M 		= 0.;
	X 		<< 0., 0., 0.;
	Xb 		<< 0., 0., 0.;
	Removed = false;
}

inline RWM_PARTICLE::RWM_PARTICLE(int type, const Vector3d& x, double m)
{
    Type	= type;
	ID		= 0;
	Tag		= 0;
	M 		= m;
	X 		= x;
	Xb 		= x;
	Removed = false;
}