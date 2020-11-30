struct SUBDOMAIN
{
	size_t Nx;		// subdomain size
	size_t Ny;
	size_t Nz;
	size_t Nc;		// total number of cells
	size_t Ncz;		// used to calculate cell index
	size_t Ncy;

	Vector3d Xmin;		// min cell position
	Vector3d Xmax;		// max cell position
};