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

template<int SType, int D>
inline void MPM<SType, D>::UpdateParticle(MPM_PARTICLE* p0)
{
	// Reset hydro force and contact force
	p0->Fh.setZero();
	p0->Fc.setZero();
	if (!p0->FixV)
	{
		Vector3d acc = Vector3d::Zero();
		Vector3d vPic = Vector3d::Zero();
		Matrix3d gradV = Matrix3d::Zero();
		double jbarDeltaJ = 0.;
		double sdf = 0.;
		for (size_t l=0; l<p0->ShapeInfo.size(); ++l)
		{
			size_t id = p0->ShapeInfo[l].ID;
			double n = p0->ShapeInfo[l].N;
			Vector3d gn = p0->ShapeInfo[l].Gn;

			acc += n*Ln[id]->F;
			vPic += n*Ln[id]->V;
			gradV += Ln[id]->V*gn.transpose();
			jbarDeltaJ += n*Ln[id]->JbarDeltaJ;
			sdf += n*Ln[id]->SDF;
		}

		Vector3d vFlip = p0->V + acc*Dt;

		p0->V = Eta*vFlip + (1.-Eta)*vPic;
		p0->X += vPic*Dt/* + 0.5*acc*Dt*Dt*/;

		p0->L = gradV;
		p0->JbarDeltaJn = jbarDeltaJ;
	}
	else
	{
		p0->V = p0->Vf;
		// p0->X += p0->V*Dt;
	}

	for (size_t d=0; d<D; ++d)
	{
		if (p0->X(d)<MinX(d))			p0->X(d) += DomSize(d);
		else if (p0->X(d)>MaxX(d))		p0->X(d) -= DomSize(d);
	}
}

template<int SType, int D>
inline void MPM<SType, D>::UpdateStress(MPM_PARTICLE* p0)
{
	// calculate deformation tensor rate
	Matrix3d deltaF = Matrix3d::Identity() + p0->L*Dt;
	double bar = 1.;
	if (UseFbar)	bar = pow(p0->JbarDeltaJn/p0->JbarDeltaJ ,1./(double) D);
	Matrix3d deltaFbar = bar*deltaF;
	p0->Td *= deltaFbar;
	// update Jbar*DeltaJ
	p0->JbarDeltaJ = p0->Td.determinant()*deltaF.determinant();

	p0->UpdatePSize(UpdatePsizeMethod);
	// Update strain
	Matrix3d de = 0.5*(deltaFbar + deltaFbar.transpose()) - Matrix3d::Identity();
	// Update stress
	Matrix3d w = 0.5*(deltaFbar - deltaFbar.transpose());
	p0->Stress += w*p0->Stress+p0->Stress*w.transpose();
	// p0->Vol 	*= 1.+de.trace();
	p0->Vol = p0->Td.determinant()*p0->Vol0;
	if (p0->MID==0)
	{
		Material::LinearElasticRateForm(de, p0->MParas[0], p0->MParas[1], p0->Stress);
	}
	else if (p0->MID==1)
	{
		double rho0 = p0->M/p0->Vol0;
		// double rho = p0->M/p0->Vol;
		double rho = rho0/p0->Td.determinant();
		p0->P = Material::EOSMonaghan(rho, rho0, Gamma, Cs);
		Material::Newtonian(de, p0->MParas[0], p0->P, p0->Stress);
	}
	else if (p0->MID==2)
	{
		Material::MohrCoulomb(de, p0->MParas[0], p0->MParas[1], p0->MParas[2], p0->MParas[3], p0->MParas[4], p0->Stress);
	}
	else if (p0->MID==3)
	{
		Material::DruckerPrager(de, p0->MParas[0], p0->MParas[1], p0->MParas[2], p0->MParas[3], p0->MParas[4], p0->MParas[5], p0->Stress);
	}
}

template<int SType, int D>
inline void MPM<SType, D>::NodeToParticle()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		UpdateParticle(Lp[p]);
		UpdateStress(Lp[p]);
	}
}