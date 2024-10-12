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

inline void MPM_PARTICLE::SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c)
{
	MID = 3;
	double mu = 0.5*young/(1.+poisson);
	double la = young*poisson/(1.+poisson)/(1.-2.*poisson);
	double A_dp = 0.;
	double B_dp = 0.;
	double Ad_dp = 0.;

	// inner
	if (dptype==0)
	{
		double bot = sqrt(3)*(3.+sin(phi));
		A_dp = 6.*sin(phi)/bot;
		B_dp = 6.*cos(phi)/bot;
		Ad_dp = 6.*sin(psi)/(sqrt(3)*(3.+sin(psi)));
	}
	// outer
	else if (dptype==1)
	{
		double bot = sqrt(3)*(3.-sin(phi));
		A_dp = 6.*sin(phi)/bot;
		B_dp = 6.*cos(phi)/bot;
		Ad_dp = 6.*sin(psi)/(sqrt(3)*(3.-sin(psi)));
	}
	// plane strain
	else if (dptype==2)
	{
		double bot = sqrt(9.+12*tan(phi)*tan(phi));
		A_dp = 3.*tan(phi)/bot;
		B_dp = 3./bot;
		Ad_dp = 3.*tan(psi)/sqrt(9.+12*tan(psi)*tan(psi));
	}
	else
	{
		cout << "not supported DP type" << endl;
		abort();
	}
	MParas[0] = mu;
	MParas[1] = la;
	MParas[2] = A_dp;
	MParas[3] = B_dp;
	MParas[4] = Ad_dp;
	MParas[5] = c;
}