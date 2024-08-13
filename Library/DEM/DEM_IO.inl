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

#ifndef DEM_IO_H
#define DEM_IO_H

inline void DEM::WriteFileH5(int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(9) << setfill('0') << n;			
	string file_name_h5 = "DEM_"+out.str()+".h5";

    hid_t     file_id;
    file_id = H5Fcreate(file_name_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);	
	
	size_t npar = Lp.size();

	hsize_t	dims_scalar[1] = {npar};			//create data space.
	hsize_t	dims_vector[1] = {3*npar};			//create data space.

	hsize_t n_points = 0;						// total number of points
	hsize_t n_faces  = 0;						// total number of faces
	hsize_t n_fe  = 0;							// total number of face elements

	for (size_t i=0; i<Lp.size(); ++i)
	{
		n_points += Lp[i]->P0.size();
		n_faces  += Lp[i]->Faces.size();
		n_fe 	 += Lp[i]->Nfe;
	}

	hsize_t	dims_points [1] = {3*n_points};			//create data space.
	hsize_t	dims_faces  [1] = {n_fe};				//create data space.
	hsize_t	dims_fscalar[1] = {n_faces};			//create data space.

	float* r_h5 	= new float[  npar];
	float* rho_h5 	= new float[  npar];
	float* vol_h5 	= new float[  npar];
	float* tag_h5 	= new float[  npar];
	float* pos_h5 	= new float[3*npar];
	float* vel_h5 	= new float[3*npar];
	float* agv_h5	= new float[3*npar];
	float* fh_h5 	= new float[3*npar];
	float* poi_h5	= new float[3*n_points];
	int* 	fac_h5	= new int   [n_fe];
	float* fv_h5 	= new float[n_faces];
	float* ftag_h5 = new float[n_faces];
	size_t count_p = 0;
	size_t count_f = 0;
	size_t count_fv = 0;

	for (size_t i=0; i<npar; ++i)
	{
		Vector3d agv = Lp[i]->Qf._transformVector(Lp[i]->W);
        r_h5  [  i  ] 	= Lp[i]->R;
        rho_h5[  i  ] 	= Lp[i]->Rho;
        vol_h5[  i  ] 	= Lp[i]->Vol;
        tag_h5[  i  ] 	= Lp[i]->Tag;
        // tag_h5[  i  ] 	= Lp[i]->ID;

		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);
		fh_h5[3*i  ] 	= Lp[i]->Fh(0);
		fh_h5[3*i+1] 	= Lp[i]->Fh(1);
		fh_h5[3*i+2] 	= Lp[i]->Fh(2);
		agv_h5[3*i  ] 	= agv(0);
		agv_h5[3*i+1] 	= agv(1);
		agv_h5[3*i+2] 	= agv(2);

		size_t add_poi = count_p;
		for (size_t j=0; j<Lp[i]->P0.size(); ++j)
		{
			Vector3d pointJ = Lp[i]->GetP(j);
			poi_h5[3*count_p  ] = pointJ(0);
			poi_h5[3*count_p+1] = pointJ(1);
			poi_h5[3*count_p+2] = pointJ(2);
			count_p++;
		}
		for (size_t k=0; k<Lp[i]->Faces.size(); ++k)
		{
			size_t num = Lp[i]->Faces[k].size()+1;
			fac_h5[count_f  ] = num;
			for (size_t m=1; m<num; ++m)
			{
				fac_h5[count_f+m] = Lp[i]->Faces[k](m-1)+add_poi;
			}

			count_f += num;
			fv_h5[count_fv] 	= Lp[i]->V.norm();
			ftag_h5[count_fv] 	= Lp[i]->Tag;
			// ftag_h5[count_fv] 	= Lp[i]->FconvID[count_fv];
			count_fv++;
		}
	}

	H5LTmake_dataset_float(file_id,"Radius",1,dims_scalar,r_h5);
	H5LTmake_dataset_float(file_id,"Rho",1,dims_scalar,rho_h5);
	H5LTmake_dataset_float(file_id,"Vol",1,dims_scalar,vol_h5);
	H5LTmake_dataset_float(file_id,"Tag",1,dims_scalar,tag_h5);
	H5LTmake_dataset_float(file_id,"Position",1,dims_vector,pos_h5);
	H5LTmake_dataset_float(file_id,"Velocity",1,dims_vector,vel_h5);
	H5LTmake_dataset_float(file_id,"AngularVelocity",1,dims_vector,agv_h5);
	H5LTmake_dataset_float(file_id,"HydroForce",1,dims_vector,fh_h5);
	H5LTmake_dataset_float(file_id,"Points",1,dims_points,poi_h5);
	H5LTmake_dataset_int(file_id,"Faces",1,dims_faces,fac_h5);
	H5LTmake_dataset_float(file_id,"FaceVelocity",1,dims_fscalar,fv_h5);
	H5LTmake_dataset_float(file_id,"FaceTag",1,dims_fscalar,ftag_h5);

	delete[] r_h5;
	delete[] rho_h5;
	delete[] tag_h5;
	delete[] vol_h5;
	delete[] pos_h5;
	delete[] vel_h5;
	delete[] agv_h5;
	delete[] fh_h5;
	delete[] poi_h5;
	delete[] fac_h5;
	delete[] fv_h5;
	delete[] ftag_h5;

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	string file_name_xmf = "DEM_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";

    if (n_faces>0)
    {
	    oss << "   <Grid Name=\"DEM_FACES\">\n";
	    oss << "     <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << n_faces << "\">\n";
	    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << n_fe << "\">\n";
	    oss << "        " << file_name_h5 <<":/Faces \n";
	    oss << "       </DataItem>\n";
	    oss << "     </Topology>\n";
	    oss << "     <Geometry GeometryType=\"XYZ\">\n";
	    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_points << " 3\" >\n";
	    oss << "        " << file_name_h5 <<":/Points \n";
	    oss << "       </DataItem>\n";
	    oss << "     </Geometry>\n";
	    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
	    oss << "       <DataItem Dimensions=\"" << n_faces << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
	    oss << "        " << file_name_h5 <<":/FaceVelocity\n";
	    oss << "       </DataItem>\n";
	    oss << "     </Attribute>\n";
	    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
	    oss << "       <DataItem Dimensions=\"" << n_faces << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
	    oss << "        " << file_name_h5 <<":/FaceTag\n";
	    oss << "       </DataItem>\n";
	    oss << "     </Attribute>\n";
	    oss << "   </Grid>\n";
    }

    oss << "   <Grid Name=\"DEM_CENTER\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << npar << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << npar << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Rho\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Rho \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Vol\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Vol \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngularVelocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/AngularVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"HydroForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << npar << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/HydroForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}

// inline void DEM::WriteContactForceFileH5(int n)
// {
// 	stringstream	out;							//convert int to string for file name.
// 	out << setw(9) << setfill('0') << n;			
// 	string file_name_h5 = "DEM_Force_"+out.str()+".h5";

// 	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
// 	hsize_t	dims_scalar[1] = {Lc.size()};			//create data space.
// 	hsize_t	dims_vector[1] = {3*Lc.size()};			//create data space.

// 	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
// 	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

// 	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
// 	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);

// 	double* p0_h5 	= new double[  Lc.size()];
// 	double* p1_h5 	= new double[  Lc.size()];
// 	double* fc_h5 	= new double[3*Lc.size()];

// 	for (size_t i=0; i<Lc.size(); ++i)
// 	{
// 		size_t p = Lc[i].PID;
// 		size_t q = Lc[i].QID;
//         p0_h5[  i  ] 	= p;
//         p1_h5[  i  ] 	= q;
// 		fc_h5[3*i  ] 	= Lp[q]->Fc(0);
// 		fc_h5[3*i+1] 	= Lp[q]->Fc(1);
// 		fc_h5[3*i+2] 	= Lp[q]->Fc(2);
// 	}

// 	DataSet	*dataset_p0 = new DataSet(file.createDataSet("P0", PredType::NATIVE_DOUBLE, *space_scalar));
// 	DataSet	*dataset_p1	= new DataSet(file.createDataSet("P1", PredType::NATIVE_DOUBLE, *space_scalar));
//     DataSet	*dataset_fc	= new DataSet(file.createDataSet("ContactForce", PredType::NATIVE_DOUBLE, *space_vector));

// 	dataset_p0->write(p0_h5, PredType::NATIVE_DOUBLE);
// 	dataset_p1->write(p1_h5, PredType::NATIVE_DOUBLE);
// 	dataset_fc->write(fc_h5, PredType::NATIVE_DOUBLE);

// 	delete space_scalar;
// 	delete space_vector;

// 	delete dataset_p0;
// 	delete dataset_p1;
// 	delete dataset_fc;

// 	delete[] p0_h5;
// 	delete[] p1_h5;
// 	delete[] fc_h5;

// 	file.close();
// }

#endif