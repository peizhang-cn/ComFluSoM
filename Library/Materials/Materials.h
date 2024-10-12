/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#pragma once

namespace Material
{
	template<typename T>
	void LinearElastic(T& e, double mu, double la, T& stress);
	template<typename T>
	void LinearElasticRateForm(T& de, double mu, double la, T& stress);
	template<typename T>
	void VKElastic(T& td, double mu, double la, T& stress);
	template<typename T>
	void YeohModel(T& td, double C10, double C20, double C30, double K, T& stress);
	template<typename T>
	void DemirayModel(T& td, double D1, double D2, double K, T& stress);
	template<typename T>
	void Newtonian(T& de, double Mu, double P, T& stress);
	void MohrCoulomb(Matrix3d& de, double Mu, double La, double Phi, double Psi, double C, Matrix3d& stress);
	void DruckerPrager(Matrix3d& de, double Mu, double La, double A_dp, double B_dp, double Ad_dp, double C, Matrix3d& stress);
	double EOSMorris(double rho, double Cs);
	double EOSMonaghan(double rho, double rho0, double Gamma, double Cs);
}

#include "LinearElastic.inl"
#include "HyperElastic.inl"
#include "Newtonian.inl"
#include "MohrCoulomb.inl"
#include "DruckerPrager.inl"
#include "EquationOfState.inl"