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

#ifndef DEM_CONTACT_H
#define DEM_CONTACT_H

inline void DEM::PeriodiParticlePosition2P(DEM_PARTICLE* pi, DEM_PARTICLE* pj, bool forP, int peri[3], Vector3d& Xi, Vector3d& Xj)
{
	Xi = pi->X;
	Xj = pj->X;

	// move particle position if periodic boundary applies
	Vector3d xperi (peri[0]*Lx, peri[1]*Ly, peri[2]*Lz);
	if (forP)	Xi += xperi;
	else		Xj += xperi;
}

inline void DEM::ContactInformation2P(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, vector<CONTACT_INFO>& lci)
{
	bool mayContact = false;

	Vector3d cpi, cpj;								// closest point on i and j, contact point
	cpi = Xi; cpj = Xj;								// set to sphere center for sphere collisions

	Vector3d n = cpi-cpj;							// Normal direction (pj pinnts to pi)
	double delta = pi->R+pj->R-n.norm(); 			// Overlapping distance for sphere

	if (delta+Hn>0.)	mayContact = true;
	else				mayContact = false;

	if (mayContact)			// when surrounding sphere is in contact
	{
		lci.clear();
		// sphere-sphere
		if ((pi->ShapeType==1 && pj->ShapeType==1))
		{
			n.normalize();									// Normalize contact normal
			CONTACT_INFO ci;
			ci.Delta = delta;
			ci.Xc = cpj+(pj->R-0.5*delta)*n;			// Contact point
			ci.Nc = n;
			lci.push_back(ci);
		}
		// sphere-cuboid
		else if ((pi->ShapeType==1 && pj->ShapeType==2) || (pi->ShapeType==2 && pj->ShapeType==1))
		{
			Sphere2Cuboid(Xi, Xj, pi, pj, lci);
		}
		// sphere-cylinder
		else if ((pi->ShapeType==1 && pj->ShapeType==3) || (pi->ShapeType==3 && pj->ShapeType==1))
		{
			Sphere2Cylinder(Xi, Xj, pi, pj, lci);
		}
		else
		{
			cout << "undefined contact shape pair" << endl;
			abort();
		}
		// // nonconvex metaball-metaball
		// else if (pi->ShapeType==4 && pj->ShapeType==4)
		// {
		// 	// cout << "in nonconvex metaball-metaball" << endl;
		// 	// if (pi->isConvex*pj->isConvex==false)	NonConvexMetaball2Metaball(Xi, Xj, pi, pj, lci);
		// }
		// // nonconvex metaball-cubiod
		// else if ((pi->ShapeType==4 && pj->ShapeType==2) || (pi->ShapeType==2 && pj->ShapeType==4))
		// {
		// 	if (pi->isConvex*pj->isConvex==false)	NonConvexMetaball2Cubiod(Xi, Xj, pi, pj, lci);
		// }
		// // nonconvex metaball-polyhedra
		// else if ((pi->ShapeType==2 && pj->ShapeType==4) || (pi->ShapeType==4 && pj->ShapeType==2))
		// {
		// 	if (pi->convex*pj->convex==false)	Metaball2PolyhedraNonConvex(Xi, Xj, pi, pj, lci);
		// }
		// // convex polyhedra-polyhedra
		// else if (pi->ShapeType==2 && pj->ShapeType==2)
		// {
		// 	if (pi->convex*pj->convex==false)
		// 	{
		// 		cout << "\033[1;31mError: only convex polyhedra-polyhedra collision is supported for now!\033[0m\n";		
		// 		exit(0);
		// 	}
		// 	else Polyhedra2PolyhedraConvex(Xi, Xj, pi, pj, lci);
		// }
		// // convex sphere-polyhedra
		// else if ((pi->ShapeType==1 && pj->ShapeType==2) || (pi->ShapeType==2 && pj->ShapeType==1))
		// {
		// 	if (pi->convex*pj->convex==false)
		// 	{
		// 		cout << "\033[1;31mError: only convex sphere-polyhedra collision is supported for now!\033[0m\n";		
		// 		exit(0);				
		// 	}
		// 	else Sphere2PolyhedraConvex(Xi, Xj, pi, pj, lci);
		// }
	}
}

inline void DEM::FrictionTangentForce(Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, double kt, double gt, Vector3d& cp, Vector3d& n, Vector3d& fn, Vector3d& ft)
{
	double mius = FsTable[pi->Material][pj->Material];
	if (mius!=0.)
	{
		// Relative velocity at the contact point
		Vector3d wi = pi->Qf._transformVector(pi->W);
		Vector3d wj = pj->Qf._transformVector(pj->W);

		Vector3d vij = pi->V + wi.cross(cp-Xi) - pj->V - wj.cross(cp-Xj);
		// Relative tangential velocity at the contact point
		Vector3d vt = vij-n.dot(vij)*n;						// eq.9
		// Update tangential spring
		Vector3d xi0 = vt*Dt;

		Vector3d xc = pi->Qfi._transformVector(cp-Xi);	// contact point under pi's body frame

		size_t key = pj->ID;
		bool convex = pi->isConvex&&pj->isConvex;
		if (convex)
		{
			if (pi->CFMap.count(key)>0)	xi0 += pi->CFMap[key];
		}
		else
		{
			if (pi->NFMap.count(key)>0)
			{
				// cout << "pi->NFMap[key].size(): " << pi->NFMap[key].size() << endl;
				size_t cid = 0;
				double cdis = (xc-get<1>(pi->NFMap[key][0])).norm();
				for (size_t i=1; i<pi->NFMap[key].size(); ++i)
				{
					double cdisi = (xc-get<1>(pi->NFMap[key][i])).norm();
					if (cdisi<cdis)
					{
						cid = i;
						cdis = cdisi;
					}
				}
				double maxDis = get<0>(pi->NFMap[key][cid]);
				// maxDis = 1.e-5;
				if (cdis<=maxDis)
				{
					xi0 += get<2>(pi->NFMap[key][cid]);
					// cout << "find tangent spring" << endl;
				}
			}
		}

		// Project to current tangential plane
		Vector3d xi = xi0 - n.dot(xi0)*n;					// eq.17
		// Static tangential force
		double fts = mius*fn.norm();
		ft = -kt*xi-gt*vt;									// eq.18
		if (ft.norm()>fts)
		{
			Vector3d t = ft.normalized();
			ft = fts*t;										// eq.21
			xi = -(fts*t+gt*vt)/kt;							// eq.20
		}

		if (convex)	pi->CFMapt[key] = xi;
		else
		{
			double maxDis = 1.1*vt.norm()*Dt;
			FRICTION_INFO fi = make_tuple(maxDis, xc, xi);
			pi->NFMapt[key].push_back(fi);
		}
	}
}

// inline void DEM::RollingResistance(DEM_PARTICLE* pi, DEM_PARTICLE* pj, double delta, double kr, double gr, Vector3d& n, Vector3d& fn, Vector3d& xir, Vector3d& armr)
// {
// 	// Reduced radius
// 	double rij = pi->R*pj->R/(pi->R+pj->R);
// 	// Rolling velocity
// 	Vector3d vij = rij*(n.cross(pj->W)-n.cross(pi->W));	// eq.15
// 	// Relative tangential velocity at the contact point
// 	Vector3d vt = vij-n.dot(vij)*n;						// eq.9
// 	// Update tangential spring
// 	Vector3d xi0 = vt*Dt;
// 	size_t key = Key(pi->ID, pj->ID);
// 	if (RMap.count(key)>0)	xi0 += RMap[key];			// eq.19
// 	// Project to current tangential plane
// 	xir = xi0 - n.dot(xi0)*n;							// eq.17
// 	// Static tangential force
// 	double fts = RsTable[pi->Material][pj->Material]*fn.norm();
// 	// Tangential force
// 	Vector3d ft = -kr*xir-gr*vt;						// eq.18
// 	if (ft.norm()>fts)
// 	{
// 		Vector3d t = ft.normalized();
// 		ft = fts*t;										// eq.21
// 		xir = -(fts*t+gr*vt)/kr;						// eq.20
// 	}
// 	// Torque with normalized arm
// 	armr = -n.cross(ft);
// }


inline void DEM::ContactForce2P(int cmType, Vector3d Xi, Vector3d Xj, DEM_PARTICLE* pi, DEM_PARTICLE* pj, vector<CONTACT_INFO>& lci)
{
	size_t lciSize = lci.size();
	if (lciSize>0)
	{
		Vector3d fc (0.,0.,0.);
		Vector3d tci (0.,0.,0.);
		Vector3d tcj (0.,0.,0.);

		for (size_t c=0; c<lci.size(); ++c)
		{
			double delta = lci[c].Delta;
			Vector3d xc = lci[c].Xc;
			Vector3d nc = lci[c].Nc;

			// get contact paramaters
			CONTACT_PARA para;
			ContactPara(cmType, pi, pj, delta, para);

			Vector3d vn = (pj->V-pi->V).dot(nc)*nc;					// Relative velocity in normal direction
			Vector3d fn = para.Kn*delta*nc + para.Gn*vn;			// Normal contact force
			// Vector3d fn = para.Kn*delta*nc;
			// cout << "delta: " << delta << endl;
			// cout << "fn: " << fn.transpose() << endl;
			// cout << "nc: " << nc.transpose() << endl;
			// cout << "vn: " << vn.transpose() << endl;
			// cout << "para.Kn*delta*nc: " << (para.Kn*delta*nc).transpose() << endl;
			// cout << "para.Gn*vn: " << (para.Gn*vn).transpose() << endl;
			// cout << "para.Kn: " << para.Kn << endl;
			// cout << "para.Gn: " << para.Gn << endl;
			// abort();
			Vector3d ft (0.,0.,0.);
			FrictionTangentForce(Xi, Xj, pi, pj, para.Kt, para.Gt, xc, nc, fn, ft);

			Vector3d fnt = fn + ft;
			fc += fnt;
			tci += pi->Qfi._transformVector((xc-Xi).cross(fnt));
			tcj -= pj->Qfi._transformVector((xc-Xj).cross(fnt));
		}
		lci.clear();

		// cout << "tcj: " << tcj.transpose() << endl;
		// cout << "I: " << pj->I.transpose() << endl;

		for (size_t d=0; d<D; ++d)
		{
			#pragma omp atomic
			pi->Fc(d) += fc(d);
			#pragma omp atomic
			pj->Fc(d) -= fc(d);
			#pragma omp atomic
			pi->Tc(d) += tci(d);
			#pragma omp atomic
			pj->Tc(d) += tcj(d);
		}
	}
}

inline void DEM::Contact2P(int cmType, DEM_PARTICLE* pi, DEM_PARTICLE* pj, bool forP, int peri[3])
{
	Vector3d Xi, Xj;
	PeriodiParticlePosition2P(pi, pj, forP, peri, Xi, Xj);
	vector<CONTACT_INFO> lci(0);
	ContactInformation2P(Xi, Xj, pi, pj, lci);
	ContactForce2P(cmType, Xi, Xj, pi, pj, lci);
}

inline void DEM::Contact(bool writeFc, int n, bool show)
{
	#pragma omp parallel for schedule(static, 1) num_threads(Nproc)
	for (size_t p=0; p<Lp.size(); ++p)
	{
		for (size_t i=0; i<Lp[p]->Lc.size(); ++i)
		{
			size_t q = Lp[p]->Lc[i].first;
			Contact2P(CMTypeIndex, Lp[p], Lp[q], false, Lp[p]->Lc[i].second);
		}
    }
}

#endif