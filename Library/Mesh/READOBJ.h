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

// read OBJ file data: vertices and faces and face normal
void ReadOBJ(string fname, vector<Vector3d>& ver, vector<VectorXi>& face, vector<Vector3d>& fnorm)
{
	cout << "Start reading from file: " << fname << endl; 
	string line;
	char c, c1;
	double x, y, z;

	ifstream in(fname);
	if (!in)
	{
		std::cout << "Could not open file: " << fname << std::endl;
		abort();
	}
	while ( getline( in, line ) )                          
	{
		istringstream ss( line );
	  	if (line[0]=='v' && line[1]==' ')
	  	{
	  		ss >> c >> x >> y >> z;
	        Vector3d v (x,y,z);
	        ver.push_back(v);
	        // cout << "v: " << v.transpose() << endl;
	  	}
	  	else if (line[0]=='v' && line[1]=='n')
	  	{
	  		ss >> c >> c1 >> x >> y >> z;
	  		Vector3d fn (x,y,z);
	  		fnorm.push_back(fn);
	  		// cout << "fn: " << fn.transpose() << endl;
	  	}
	  	else if (line[0]=='f')
	  	{
	  		vector<string> ls;
			string lse;
			while (ss>>lse)	ls.push_back(lse);
	        VectorXi f;
	        f.resize(ls.size()-1);
	        for (size_t m=1; m<ls.size(); ++m)
	        {
	        	f(m-1) = stoi(ls[m])-1;
	        }
	        face.push_back(f);
	        // cout << "f: " << f.transpose() << endl;
	  	}
	}
	in.close();
	cout << "Reading mesh finshed." << endl;
	// abort();
}