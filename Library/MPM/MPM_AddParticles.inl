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
inline void MPM<SType, D>::AddParticle(int tag, Vector3d& x, double m)
{
    Lp.push_back(new MPM_PARTICLE(ReservedN, tag,x,m));
    Lp[Lp.size()-1]->ID = Lp.size()-1;
}

template<int SType, int D>
inline void MPM<SType, D>::AddBoxParticles(int tag, Vector3d& x0, Vector3d& x1, double dpx, double rho)
{
	Vector3d l = x1-x0;
	Vector3i maxx = Vector3i::Zero();

	for (size_t d=0; d<D; ++d)
	{
		maxx(d) = round(l(d)/dpx);
		if (abs(l(d)-maxx(d)*dpx)>(1.e-12*dpx))
		{
			cout << l(d) << endl;
			cout << maxx(d)*dpx << endl;
			cout << abs(l(d)-maxx(d)*dpx) << endl;
			cout << "\033[1;31mError: When devide length of the box by dpx, remainder is not zero! \033[0m\n";		
			exit(0);
		}
		maxx(d) -= 1;
	}

	double vol = pow(dpx,D);
	double mp = vol*rho;

	for (int k=0; k<=maxx(2); ++k)
	for (int j=0; j<=maxx(1); ++j)
	for (int i=0; i<=maxx(0); ++i)
	{
		Vector3d x = Vector3d::Zero();
						x(0) = dpx*(i + 0.5)+x0(0);
		if (D>1)		x(1) = dpx*(j + 0.5)+x0(1);
		if (D>2)		x(2) = dpx*(k + 0.5)+x0(2);

		AddParticle(tag, x, mp);
		size_t p = Lp.size()-1;
		if (SType==1)
		{
			for (size_t d=0; d<D; ++d)
			{
				Lp[p]->PSize0(d) = 0.5*dpx;
				Lp[p]->PSize(d) = Lp[p]->PSize0(d);
			}
		}
    	Lp[p]->Vol0 = vol;
    	Lp[p]->Vol 	= vol;
	}
}

template<int SType, int D>
inline void MPM<SType, D>::DeleteParticles()
{
	vector <MPM_PARTICLE*>	Lpt;
	Lpt.resize(0);

	vector <MPM_PARTICLE*>	Lrm;
	Lrm.resize(0);

	for (size_t p=0; p<Lp.size(); ++p)
	{
		if (!Lp[p]->Removed)	Lpt.push_back(Lp[p]);
		else					Lrm.push_back(Lp[p]);
	}
	Lp = Lpt;

	for (size_t p=0; p<Lp.size(); ++p)
	{
		Lp[p]->ID = p;
	}

	// for (size_t p=0; p<Lrm.size(); ++p)
	// {
	// 	delete[] Lrm[p];
	// }

}