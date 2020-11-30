#include "../HEADER.h"
#include <MPM.h>

class MMPM
{
public:
	MMPM();
	~MMPM();
	MMPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx);
	Init();

	MPM*				DomF;				// Domain of fluid
	MPM*				DomS;				// Domain of Soil

	double 				Rhosp;				// density of soil particles
};

void MMPM::MMPM(size_t ntype, size_t nx, size_t ny, size_t nz, Vector3d dx)
{
	DomF = MPM(ntype, nx, ny, nz, dx);
	DomS = MPM(ntype, nx, ny, nz, dx);
}

void MMPM::Init()
{
	DomF->Init();
	DomS->Init();
}

void MMPM::CalPorosity()
{
	double poros = 1.-DomS->Ln[id]->M/Rhosp;
}

void MMPM::SolveMUSL(int tt, int ts)
{
	DomF->ParticleToNode();
	DomS->ParticleToNode();
}