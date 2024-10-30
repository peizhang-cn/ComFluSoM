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
#pragma once

// calculate mass center, moment of inertia and other properties of grouped spheres/disks
namespace GroupedSpheresProperties
{
	template<typename T>
	void CalProperties2D(vector<T> pos, vector<T> r, double rho, double& vol, double& m, T& xc, double& inertia)
	{
		vol = 0.;
		m = 0.;
		xc.setZero();
		inertia = 0.;
		for (size_t i=0; i<pos.size(); ++i)
		{
			double voli = M_PI*r[i]*r[i];
			vol += voli;
			double mi = rho*voli;
			m += mi;
			xc += pos[i]*voli;
			// 这里的参考点是pos[0]，最后需要都移动到质心
			inertia += 0.5*mi*(r[i]*r[i]) + mi*(pos[i]-pos[0]).squaredNorm();
		}

		// correct for overlap
		for (size_t i=0; i<pos.size(); ++i)
		for (size_t j=i+1; j<pos.size(); ++j)
		{
			double dis = r[i]+r[j]-(pos[i]-pos[j]).norm();
			if (dis>0.)
			{
				// 把重叠部分分为两个cap,分别计算
				// 需要两种情况，一种小球球心不在大球内，一种小球在大球内
				// 小球在大球内时，可以计算小球较小部分的质心、惯性矩，然后反推较大部分的
				// 修正质心、体积、质量和惯性矩
			}
		}

		xc /= vol;
	}

	template<typename T>
	void MonteCarlo2D(size_t n, vector<T>& pos, vector<double>& r, double rho, double& area, double& m, T& xc, double& inertia)
	{
		double minx = numeric_limits<double>::max();
		double maxx = -numeric_limits<double>::max();
		double miny = numeric_limits<double>::max();
		double maxy = -numeric_limits<double>::max();

		for (size_t i=0; i<pos.size(); ++i)
		{
			minx = min(minx, pos[i][0]-r[i]);
			maxx = max(maxx, pos[i][0]+r[i]);
			miny = min(miny, pos[i][1]-r[i]);
			maxy = max(maxy, pos[i][1]+r[i]);
		}
		double lx = maxx - minx;
		double ly = maxy - miny;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0., 1.);

		double sum_x = 0.0, sum_y = 0.0;
		double sum_r2 = 0.0;

		size_t count = 0;
		for (size_t i=0; i<n; ++i)
		{
			double xi = minx + lx*dis(gen);
			double yi = miny + ly*dis(gen);

			bool isInside = false;
			for (size_t j=0; j<pos.size(); ++j)
			{
				double dis = r[j]*r[j] - pow(pos[j][0]-xi,2) - pow(pos[j][1]-yi,2);
				if (dis>0.)
				{
					isInside = true;
					count++;
					break;
				}
			}
			if (isInside)
			{
				sum_x += xi;
				sum_y += yi;
				sum_r2 += pow(xi-minx,2) + pow(yi-miny,2);
			}
		}

		area = count*(lx*ly)/n;
		m = rho*area;
		xc[0] = sum_x/count;
		xc[1] = sum_y/count;

		inertia = m/count*sum_r2;
		inertia -= m*(pow(xc[0]-minx, 2) + pow(xc[1]-miny, 2));
	}

	void MonteCarlo3D(size_t n, vector<Vector3d>& pos, vector<double>& r, double rho, double& vol, double& m, Vector3d& xc, Matrix3d& inertia)
	{
		Vector3d min0, max0;
		for (size_t d=0; d<3; ++d){
			min0[d] = numeric_limits<double>::max();
			max0[d] = -numeric_limits<double>::max();
		}

		for (size_t i=0; i<pos.size(); ++i){
			for (size_t d=0; d<3; ++d){
				min0[d] = min(min0[i], pos[i][d]-r[i]);
				max0[d] = max(max0[i], pos[i][d]+r[i]);
			}
		}
		Vector3d l0 = max0 - min0;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0., 1.);

		Vector3d xc0 (0., 0., 0.);
		Matrix3d inertia0 = Matrix3d::Zero(); 

		size_t count = 0;
		for (size_t i=0; i<n; ++i)
		{
			Vector3d xi;
			for (size_t d=0; d<3; ++d){
				xi[d] = min0[d] + l0[d]*dis(gen);
			}

			bool isInside = false;
			for (size_t j=0; j<pos.size(); ++j)
			{
				double dis = r[j]*r[j] - (pos[j]-xi).squaredNorm();
				if (dis>0.)
				{
					isInside = true;
					count++;
					break;
				}
			}
			if (isInside)
			{
				xc0 += xi;
				xi -= min0;

				inertia0(0, 0) += xi(1) * xi(1) + xi(2) * xi(2);
				inertia0(1, 1) += xi(0) * xi(0) + xi(2) * xi(2);
				inertia0(2, 2) += xi(0) * xi(0) + xi(1) * xi(1);

				inertia0(0, 1) -= xi(0) * xi(1);
				inertia0(0, 2) -= xi(0) * xi(2);
				inertia0(1, 2) -= xi(1) * xi(2);

				inertia0(1, 0) -= xi(1) * xi(0); // Same as (0, 1), but it's fine to do both
				inertia0(2, 0) -= xi(2) * xi(0); // Same as (0, 2)
				inertia0(2, 1) -= xi(2) * xi(1); // Same as (1, 2)
			}
		}

		vol = count*l0[0]*l0[1]*l0[2]/n;
		m = rho*vol;
		xc = xc0/count;
		inertia0 *= m/count;

		Vector3d xcd = min0-xc;

		inertia0(0, 0) -= m * (xcd(1) * xcd(1) + xcd(2) * xcd(2));
		inertia0(1, 1) -= m * (xcd(0) * xcd(0) + xcd(2) * xcd(2));
		inertia0(2, 2) -= m * (xcd(0) * xcd(0) + xcd(1) * xcd(1));

		inertia0(0, 1) += m * xcd(0) * xcd(1);
		inertia0(0, 2) += m * xcd(0) * xcd(2);
		inertia0(1, 2) += m * xcd(1) * xcd(2);

		inertia0(1, 0) += m * xcd(1) * xcd(0); // Same as (0, 1)
		inertia0(2, 0) += m * xcd(2) * xcd(0); // Same as (0, 2)
		inertia0(2, 1) += m * xcd(2) * xcd(1); // Same as (1, 2)

		inertia = inertia0;
	}

	template<typename T>
	void Grid2D(double dx, vector<T>& pos, vector<double>& r, double rho, double& area, double& m, T& xc, double& inertia)
	{
		T min0 = T::Zero();
		T max0 = T::Zero();

		for (size_t d=0; d<2; ++d){
			min0[d] = numeric_limits<double>::max();
			max0[d] = -numeric_limits<double>::max();
		}

		for (size_t i=0; i<pos.size(); ++i){
			for (size_t d=0; d<2; ++d){
				min0[d] = min(min0[d], pos[i][d]-r[i]);
				max0[d] = max(max0[d], pos[i][d]+r[i]);
			}
		}

		T l0 = max0 - min0;
		T xc0 = T::Zero();
		double inertia0 = 0.; 

		size_t count = 0;
		size_t nx = size_t(l0[0]/dx) + 1;
		size_t ny = size_t(l0[1]/dx) + 1;

		cout << "nx: " << nx << endl;
		cout << "ny: " << ny << endl;

		for (size_t i=0; i<nx; ++i)
		for (size_t j=0; j<ny; ++j)
		{
			T xi;
			xi[0] = min0[0] + i*dx;
			xi[1] = min0[1] + j*dx;

			bool isInside = false;
			for (size_t s=0; s<pos.size(); ++s)
			{
				double dis = r[s]*r[s] - (pos[s]-xi).squaredNorm();
				if (dis>0.)
				{
					isInside = true;
					count++;
					break;
				}
			}
			if (isInside)
			{
				xc0 += xi;
				xi -= min0;
				inertia0 += xi.squaredNorm();
			}
		}

		area = count*dx*dx;
		m = rho*area;
		xc = xc0/count;

		inertia0 *= m/count;
		inertia0 -= m*(xc-min0).squaredNorm();
		inertia = inertia0;
	}

	void Grid3D(double dx, vector<Vector3d>& pos, vector<double>& r, double rho, double& vol, double& m, Vector3d& xc, Matrix3d& inertia)
	{
		Vector3d min0, max0;
		for (size_t d=0; d<3; ++d){
			min0[d] = numeric_limits<double>::max();
			max0[d] = -numeric_limits<double>::max();
		}

		for (size_t i=0; i<pos.size(); ++i){
			for (size_t d=0; d<3; ++d){
				min0[d] = min(min0[d], pos[i][d]-r[i]);
				max0[d] = max(max0[d], pos[i][d]+r[i]);
			}
		}
		Vector3d l0 = max0 - min0;

		Vector3d xc0 (0., 0., 0.);
		Matrix3d inertia0 = Matrix3d::Zero(); 

		size_t count = 0;
		size_t nx = size_t(l0[0]/dx) + 1;
		size_t ny = size_t(l0[1]/dx) + 1;
		size_t nz = size_t(l0[2]/dx) + 1;

		for (size_t i=0; i<nx; ++i)
		for (size_t j=0; j<ny; ++j)
		for (size_t k=0; k<nz; ++k)
		{
			Vector3d xi = Vector3d(i*dx, j*dx, k*dx) + min0;

			bool isInside = false;
			for (size_t s=0; s<pos.size(); ++s)
			{
				double dis = r[s]*r[s] - (pos[s]-xi).squaredNorm();
				if (dis>0.)
				{
					isInside = true;
					count++;
					break;
				}
			}
			if (isInside)
			{
				xc0 += xi;
				xi -= min0;

				inertia0(0, 0) += xi(1) * xi(1) + xi(2) * xi(2);
				inertia0(1, 1) += xi(0) * xi(0) + xi(2) * xi(2);
				inertia0(2, 2) += xi(0) * xi(0) + xi(1) * xi(1);

				inertia0(0, 1) -= xi(0) * xi(1);
				inertia0(0, 2) -= xi(0) * xi(2);
				inertia0(1, 2) -= xi(1) * xi(2);

				inertia0(1, 0) -= xi(1) * xi(0); // Same as (0, 1), but it's fine to do both
				inertia0(2, 0) -= xi(2) * xi(0); // Same as (0, 2)
				inertia0(2, 1) -= xi(2) * xi(1); // Same as (1, 2)
			}
		}

		vol = count*dx*dx*dx;
		m = rho*vol;
		xc = xc0/count;
		inertia0 *= m/count;

		Vector3d xcd = min0-xc;

		inertia0(0, 0) -= m * (xcd(1) * xcd(1) + xcd(2) * xcd(2));
		inertia0(1, 1) -= m * (xcd(0) * xcd(0) + xcd(2) * xcd(2));
		inertia0(2, 2) -= m * (xcd(0) * xcd(0) + xcd(1) * xcd(1));

		inertia0(0, 1) += m * xcd(0) * xcd(1);
		inertia0(0, 2) += m * xcd(0) * xcd(2);
		inertia0(1, 2) += m * xcd(1) * xcd(2);

		inertia0(1, 0) += m * xcd(1) * xcd(0); // Same as (0, 1)
		inertia0(2, 0) += m * xcd(2) * xcd(0); // Same as (0, 2)
		inertia0(2, 1) += m * xcd(2) * xcd(1); // Same as (1, 2)

		inertia = inertia0;
	}
}