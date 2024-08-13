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

// Add a 2d disk
inline void DEM::AddDisk2D(int tag, double r, Vector3d& x, double rho)
{
	if (D!=2)
	{
		cout << "\033[1;31mError: Disk2D only works on 3D!\033[0m\n";		
		exit(0);		
	}
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = Lp.size()-1;
	Lp[p]->SetDisk2D(r);
	Bins->AddParticleToBin(Lp[p]);
    if (Lp[p]->R>0.5*Bins->MinL)
    {
        cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";     
        exit(0);    
    }
}
// Add a sphere
inline void DEM::AddSphere(int tag, double r, Vector3d& x, double rho)
{
	if (D!=3)
	{
		cout << "\033[1;31mError: Sphere only works on 3D!\033[0m\n";		
		exit(0);		
	}
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	Lp[p]->SetSphere(r);
	Bins->AddParticleToBin(Lp[p]);
    if (Lp[p]->R>0.5*Bins->MinL)
    {
        cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";     
        exit(0);    
    }
}
// Add a metaball
inline void DEM::AddMetaball(int tag, double dis, vector<Vector3d>& metaP, VectorXd& metaK, double rho)
{
	Vector3d x (0,0,0);
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	// Lp[p]->SetMetaball(rs, metaP, metaK, UnitSpherePoints, UnitSphereFaces);
	cout << "SetMetaball start" << endl;
	Lp[p]->SetMetaball(dis, metaP, metaK);
	cout << "SetMetaball finished" << endl;
	Bins->AddParticleToBin(Lp[p]);
	if (Lp[p]->R>0.5*Bins->MinL)
	{
		cout << "\033[1;31mError: The particle size is:" << Lp[p]->R << ", larger than the bin size!\033[0m\n";		
		exit(0);	
	}
}

// // Add a ametaball
// inline void DEM::AddAMetaball(int tag, double rs, vector<Vector3d>& metaP, vector<Matrix3d>& metaW, VectorXd& metaK, double rho)
// {
// 	Vector3d x (0,0,0);
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	size_t p = Lp.size()-1;
// 	Lp[p]->ID = p;
// 	Lp[p]->SetAMetaball(rs, metaP, metaW, metaK, UnitSpherePoints, UnitSphereFaces);
// 	AddParticleToBin(Lp[p], Bins);
// 	Lp[p]->R += Hn;
// 	if (Lp[p]->R>0.5*Bins->MinL)
// 	{
// 		cout << "\033[1;31mError: The particle size is:" << Lp[p]->R << ", larger than the bin size!\033[0m\n";		
// 		exit(0);	
// 	}
// }

inline void DEM::AddCuboid(int tag, Vector3d& l, Vector3d& x, double rho)
{
	// if (D!=3)
	// {
	// 	cout << "\033[1;31mError: Cuboid only works on 3D!\033[0m\n";		
	// 	exit(0);		
	// }
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	Lp[Lp.size()-1]->SetCuboid(l(0), l(1), l(2));
	if (Lp[p]->R>0.5*Bins->MinL)
	{
		// cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
		// exit(0);
		Lw.push_back(p);
		Lp[p]->isBig = true;
	}
	else
	{
		Bins->AddParticleToBin(Lp[p]);
	}
}

inline void DEM::AddCylinder(int tag, double h, double r, Vector3d& n, Vector3d& x, double rho)
{
	if (D!=3)
	{
		cout << "\033[1;31mError: Cylinder only works on 3D!\033[0m\n";		
		exit(0);		
	}
	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
	size_t p = Lp.size()-1;
	Lp[p]->ID = p;
	Lp[Lp.size()-1]->SetCylinder(h, r, n);
	if (Lp[p]->R>0.5*Bins->MinL)
	{
		cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
		// exit(0);
		Lw.push_back(p);
		Lp[p]->isBig = true;
	}
	else
	{
		Bins->AddParticleToBin(Lp[p]);
	}
}

// inline void DEM::AddTriangle2D(int tag, vector<Vector3d> ver, double rho)
// {
// 	if (D!=2)
// 	{
// 		cout << "\033[1;31mError: Triangle2D only works on 2D!\033[0m\n";		
// 		exit(0);
// 	}
// 	if (ver.size()!=3)
// 	{
// 		cout << "\033[1;31mError: Vertices number is not 3, you sure it's Triangle2D?\033[0m\n";		
// 		exit(0);		
// 	}
// 	Vector3d x (0,0,0);
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	size_t p = Lp.size()-1;
// 	Lp[p]->ID = p;
// 	cout << "SetTriangle2D" << endl;
// 	Lp[p]->SetTriangle2D(ver);
// 	cout << "SetTriangle2D Finish" << endl;
// 	AddParticleToBin(Lp[p], Bins);
// 	if (Lp[p]->R>0.5*Bins->MinL)
// 	{
// 		cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
// 		exit(0);
// 	}
// }

// inline void DEM::AddPolygon2D(int tag, vector<Vector3d> ver, double rho)
// {
// 	if (D!=2)
// 	{
// 		cout << "\033[1;31mError: Polygon2D only works on 2D!\033[0m\n";		
// 		exit(0);
// 	}
// 	Vector3d x (0,0,0);
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	size_t p = Lp.size()-1;
// 	Lp[p]->ID = p;
// 	// cout << "SetPolygon2D" << endl;
// 	Lp[p]->SetPolygon2D(ver);
// 	// cout << "SetPolygon2D Finish" << endl;
// 	AddParticleToBin(Lp[p], Bins);
// 	if (Lp[p]->R>0.5*Bins->MinL)
// 	{
// 		cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
// 		exit(0);
// 	}
// }

// inline void DEM::AddPolygon3D(int tag, double r, vector<Vector3d> ver, double rho)
// {
// 	// if (D!=3)
// 	// {
// 	// 	cout << "\033[1;31mError: Polygon3D only works on 3D!\033[0m\n";		
// 	// 	exit(0);
// 	// }
// 	Vector3d x (0,0,0);
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	size_t p = Lp.size()-1;
// 	Lp[p]->ID = p;
// 	Lp[p]->R = r;
// 	// cout << "SetPolygon3D" << endl;
// 	Lp[p]->SetPolygon3D(ver);
// 	// cout << "SetPolygon3D Finish" << endl;
// 	Lp[p]->Fix();
// 	WallsBin wall;
// 	wall.WID = p;
// 	// // find vertices under global frame
// 	// size_t npoints = Lp[p]->P0.size();
// 	// vector<Vector3d> ver;
// 	// ver,resize(npoints);
// 	// for(size_t i=0; i<npoints; ++i)	ver[i] = Lp[p]->GetP(i);
// 	Bins->FindBinIDForPolygon3D(Lp[p]->R, Lp[p]->Normal, ver, Lp[p]->Edges, wall.Lb);
// 	Lw.push_back(wall);
// }

// inline void DEM::AddTetrahedron(int tag, vector<Vector3d> ver, double rho)
// {
// 	if (D!=3)
// 	{
// 		cout << "\033[1;31mError: Tetrahedron only works on 3D!\033[0m\n";		
// 		exit(0);		
// 	}
// 	if (ver.size()!=4)
// 	{
// 		cout << "\033[1;31mError: Vertices number is not 4, you sure it's tetrahedron?\033[0m\n";		
// 		exit(0);		
// 	}
// 	Vector3d x (0,0,0);
// 	Lp.push_back(new DEM_PARTICLE(tag, x, rho));
// 	size_t p = Lp.size()-1;
// 	Lp[p]->ID = p;
// 	cout << "SetTetrahedron" << endl;
// 	Lp[p]->SetTetrahedron(ver);
// 	cout << "SetTetrahedron Finish" << endl;
// 	AddParticleToBin(Lp[p], Bins);
// 	if (Lp[p]->R>0.5*Bins->MinL)
// 	{
// 		cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
// 		exit(0);	
// 	}
// }

inline void DEM::AddNSpheres(int tag, size_t seed, size_t np, Vector3d& x0, Vector3d& x1, double r, double surDis, double rho)
{
	// seed: seed for random numbers, np: numbers of spheres, x0,x1: min and max range of box,
	// r: radius, surDis: smallest surface distance, rho: density
	if (r>0.5*Bins->MinL)
	{
		cout << "\033[1;31mError: The particle size is larger than the bin size!\033[0m\n";		
		exit(0);	
	}
    srand (seed);
    size_t count = 0;
    Vector3d l = x1-x0;					// move x0 to origin
    for (int n=0; n<1.0e200; n++)
    {
        Vector3d pos(0.,0.,0.);
        for (size_t d=0; d<D; ++d)	pos(d) = (l(d)-2.*r)*((double)rand()/RAND_MAX)+r;
        pos += x0;
        bool boolContact = false;

        size_t bid = Bins->FindBinID(pos);		// the bin id for pos
        for (size_t m=0; m<Bins->NeiV.size(); ++m)
        {
        	size_t b = bid+Bins->NeiV[m];		// neighbor bin ID
        	if (b<Bins->Nb)
        	{
	        	BIN_LC* bin = Bins->Lb[b];			// current bin
	        	size_t np = bin->Lp.size();			// numbers of particle in this bin
	        	for (size_t i=0; i<np; ++i)			// for every particle in this bin
	        	{
	        		size_t p = bin->Lp[i];			// ID of current particle
	        		double dis = (Lp[p]->X-pos).norm()-r-Lp[p]->R;
	        		if (dis<surDis)
	        		{
	        			boolContact = true;
	        			break;
	        		}
	        	}
        	}
        }
        if (boolContact==false)   
        {
            if (D==3)		AddSphere(tag, r, pos, rho);
            else if (D==2)	AddDisk2D(tag, r, pos, rho);
            size_t p = Lp.size()-1;
            Lp[p]->UpdateBox(D);
            count++;
			// cout << "count = " << count << endl;
			// cout <<"===================" << endl;
			size_t showNp = round(0.1*((double) np));
			if (count%showNp==0){
				cout<<"Add "<< 10.*(count/showNp) <<"% of " << np << " particles"<<endl;
			}
        }
		// size_t showNp = round(0.1*((double) np));
		// if (count%showNp==0){
		// 	cout << "showNp = " << showNp << endl;
		// 	cout<<"Add "<< count/showNp <<" % particles"<<endl;
		// 	cout << "count = " << count << endl;
		// 	cout <<"===================" << endl;
		// }
        if (count>=np)   break;
    }
}

inline void DEM::AddNMetaballs(int tag, size_t seed, size_t np, Vector3d& x0, Vector3d& x1, double rho, double dis, vector<Vector3d>& metaP, VectorXd& metaK)
{
    srand (seed);
    Vector3d l = x1-x0;

    for (size_t s=0; s<np; ++s)
    {
		Vector3d x (0,0,0);
		Lp.push_back(new DEM_PARTICLE(tag, x, rho));
		size_t p = Lp.size()-1;
		Lp[p]->ID = p;
		// Lp[p]->SetMetaball(rs, metaP, metaK, UnitSpherePoints, UnitSphereFaces);
		Lp[p]->SetMetaball(dis, metaP, metaK);
    	double r = Lp[p]->R;	// Radius
		if (r>0.5*Bins->MinL)
		{
			cout << "\033[1;31mError: The particle size is:" << r << ", larger than the bin size!\033[0m\n";		
			exit(0);	
		}
    	for (int n=0; n<1.0e300; ++n)
    	{
	        Vector3d pos(0.,0.,0.);
	        for (size_t d=0; d<D; ++d)	pos(d) = (l(d)-2.*r)*((double)rand()/RAND_MAX)+r;
	        pos += x0;
	        bool boolContact = false;
	        size_t bid = Bins->FindBinID(pos);		// the bin id for pos
	        for (size_t m=0; m<Bins->NeiV.size(); ++m)
	        {
	        	size_t b = bid+Bins->NeiV[m];		// neighbor bin ID
	        	if (b<Bins->Nb)
	        	{
		        	BIN_LC* bin = Bins->Lb[b];			// current bin
		        	size_t np = bin->Lp.size();			// numbers of particle in this bin
		        	for (size_t i=0; i<np; ++i)			// for every particle in this bin
		        	{
		        		size_t p = bin->Lp[i];			// ID of current particle
		        		double dis = (Lp[p]->X-pos).norm()-r-Lp[p]->R;
		        		if (dis<0.)
		        		{
		        			boolContact = true;
		        			break;
		        		}
		        	}
	        	}
	        }
			if (!boolContact)
			{
				Vector3d axe;
				axe(0) = (double)rand()/RAND_MAX;
				axe(1) = (double)rand()/RAND_MAX;
				axe(2) = (double)rand()/RAND_MAX;
				axe.normalize();
				double ang = ((double)rand()/RAND_MAX)*2.*M_PI;

				Matrix3d m;
				m = AngleAxisd(ang, axe);
				Quaterniond q0(m);
				Lp[p]->Q0 = q0;
				Lp[p]->Qf = Lp[p]->Q0*Lp[p]->Q;	// final rotation (from object frame to lab frame)
				Lp[p]->Qfi = Lp[p]->Qf.inverse();
				Lp[p]->X = pos;
				// for (size_t i=0; i<Lp[p]->MetaP.size(); ++i)
				// {
				// 	Lp[p]->MetaP[i] = q0._transformVector(Lp[p]->MetaP0[i]);
				// 	Lp[p]->MetaP[i] += pos;
				// }
				// for (size_t i=0; i<Lp[p]->P.size(); ++i)
				// {
				// 	Lp[p]->P[i] = q0._transformVector(Lp[p]->P0[i]);
				// 	Lp[p]->P[i] += pos;
				// }
				Lp[p]->UpdateBox(D);
				Bins->AddParticleToBin(Lp[p]);
				break;
			}
    	}
    }
}