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

#ifndef DEM_CONTACT_PARAMETERS_H
#define DEM_CONTACT_PARAMETERS_H

inline double DEM::EffectiveValue(double ai, double aj)
{
	double a = 0.;
	if (ai>1.0e-12 && aj>1.0e-12)	a = ai*aj/(ai+aj);
	return (a);
}

inline void DEM::LinearContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, CONTACT_PARA& para)
{
	para.Kn 	= EffectiveValue(pi->Kn, pj->Kn);
	para.Kt 	= EffectiveValue(pi->Kt, pj->Kt);
	para.Gn  	= EffectiveValue(pi->Gn, pj->Gn);
	para.Gt 	= EffectiveValue(pi->Gt, pj->Gt);
}

inline void DEM::LinearContactParaFromCr(DEM_PARTICLE* pi, DEM_PARTICLE* pj, CONTACT_PARA& para)
{
	para.Kn 	= EffectiveValue(pi->Kn, pj->Kn);
	para.Kt 	= EffectiveValue(pi->Kt, pj->Kt);

	double me 	= EffectiveValue(pi->M, pj->M);
	// // For collision with wall
	// if (pi->ShapeType==0)		me = pj->M;
	// else if (pj->ShapeType==0)	me = pi->M;

	double lgcr = log(CrTable[pi->Material][pj->Material]);
	double beta = lgcr/sqrt(lgcr*lgcr + M_PI*M_PI);
	double b = -2.*sqrt(me)*beta;
	para.Gn = b*sqrt(para.Kn);
	para.Gt = b*sqrt(para.Kt);
}

inline void DEM::HertzContactPara(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, CONTACT_PARA& para)
{
	// EDEM 2.6 Theory Reference Guide
	double re 	= EffectiveValue(pi->R, pj->R);
	double me 	= EffectiveValue(pi->M, pj->M);

	double youngi = pi->Kn;
	double youngj = pj->Kn;
	double poissoni = pi->Gn;
	double poissonj = pj->Gn;

	double ee 	= 1./((1.-poissoni*poissoni)/youngi + (1.-poissonj*poissonj)/youngj);
	double ge 	= 1./(2.*(2.-poissoni)*(1.+poissoni)/youngi + 2.*(2.-poissonj)*(1.+poissonj)/youngj);
	// // For collision with wall
	// if (pi->ShapeType==0)
	// {
	// 	me = pj->M;
	// 	re = pj->R;
	// }
	// else if (pj->ShapeType==0)
	// {
	// 	me = pi->M;
	// 	re = pi->R;		
	// }
	double rd = sqrt(re*delta);
	double sn 	= 2.*ee*rd;
	double st 	= 8.*ge*rd;
	double lgcr = log(CrTable[pi->Material][pj->Material]);
	double beta = lgcr/sqrt(lgcr*lgcr + M_PI*M_PI);
	para.Kn 	= 4./3.*ee*rd;
	para.Kt 	= 8.*ge*rd;
	double b = -2.*sqrt(5./6.*me)*beta;
	para.Gn 	= b*sqrt(sn);
	para.Gt 	= b*sqrt(st);
}

#endif