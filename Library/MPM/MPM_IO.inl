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
inline void MPM<SType, D>::WriteFileH5Particle(string name, int n)
{
	stringstream	out;							//convert int to string for file name.
	out << setw(6) << setfill('0') << n;			
	string file_name_h5 = name+"_MPM_PARTICLE_"+out.str()+".h5";

	H5File	file(file_name_h5, H5F_ACC_TRUNC);		//create a new hdf5 file.
	
	hsize_t	dims_scalar[1] = {Lp.size()};			//create data space.
	hsize_t	dims_vector[1] = {3*Lp.size()};			//create data space.
	hsize_t	dims_tensor[1] = {6*Lp.size()};
	hsize_t	dims_true_tensor[1] = {9*Lp.size()};

	int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
	int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);
	int rank_tensor = sizeof(dims_tensor) / sizeof(hsize_t);
	int rank_true_tensor = sizeof(dims_true_tensor) / sizeof(hsize_t);

	DataSpace	*space_scalar = new DataSpace(rank_scalar, dims_scalar);
	DataSpace	*space_vector = new DataSpace(rank_vector, dims_vector);
	DataSpace	*space_tensor = new DataSpace(rank_tensor, dims_tensor);
	DataSpace	*space_true_tensor = new DataSpace(rank_true_tensor, dims_true_tensor);

	float* tag_h5 	= new float[  Lp.size()];
	float* m_h5 	= new float[  Lp.size()];
	float* vol_h5 	= new float[  Lp.size()];
	float* pos_h5 	= new float[3*Lp.size()];
	float* vel_h5 	= new float[3*Lp.size()];
	float* psize_h5	= new float[3*Lp.size()];
	float* s_h5 	= new float[6*Lp.size()];
	float* td_h5 	= new float[9*Lp.size()];

	float* pressure_h5 	= new float[  Lp.size()];

	for (size_t i=0; i<Lp.size(); ++i)
	{
		tag_h5[  i  ] 	= Lp[i]->Tag;
		// tag_h5[  i  ] 	= Lp[i]->ID;
		m_h5  [  i  ] 	= Lp[i]->M;
		vol_h5[  i  ] 	= Lp[i]->Vol;
		pos_h5[3*i  ] 	= Lp[i]->X(0);
		pos_h5[3*i+1] 	= Lp[i]->X(1);
		pos_h5[3*i+2] 	= Lp[i]->X(2);
		vel_h5[3*i  ] 	= Lp[i]->V(0);
		vel_h5[3*i+1] 	= Lp[i]->V(1);
		vel_h5[3*i+2] 	= Lp[i]->V(2);
		psize_h5[3*i  ] = Lp[i]->PSize(0);
		psize_h5[3*i+1] = Lp[i]->PSize(1);
		psize_h5[3*i+2] = Lp[i]->PSize(2);

		s_h5  [6*i  ] 	= Lp[i]->Stress(0,0);
		s_h5  [6*i+1] 	= Lp[i]->Stress(0,1);
		s_h5  [6*i+2] 	= Lp[i]->Stress(0,2);
		s_h5  [6*i+3] 	= Lp[i]->Stress(1,1);
		s_h5  [6*i+4] 	= Lp[i]->Stress(1,2);
		s_h5  [6*i+5] 	= Lp[i]->Stress(2,2);

		td_h5  [9*i  ] 	= Lp[i]->Td(0,0);
		td_h5  [9*i+1] 	= Lp[i]->Td(0,1);
		td_h5  [9*i+2] 	= Lp[i]->Td(0,2);
		td_h5  [9*i+3] 	= Lp[i]->Td(1,0);
		td_h5  [9*i+4] 	= Lp[i]->Td(1,1);
		td_h5  [9*i+5] 	= Lp[i]->Td(1,2);
		td_h5  [9*i+6] 	= Lp[i]->Td(2,0);
		td_h5  [9*i+7] 	= Lp[i]->Td(2,1);
		td_h5  [9*i+8] 	= Lp[i]->Td(2,2);

		pressure_h5[i] 		= -(Lp[i]->Stress(0,0)+Lp[i]->Stress(1,1)+Lp[i]->Stress(2,2))/3.;
	}

	DataSet	*dataset_tag	= new DataSet(file.createDataSet("Tag", PredType::NATIVE_FLOAT, *space_scalar));
	DataSet	*dataset_m		= new DataSet(file.createDataSet("Mass", PredType::NATIVE_FLOAT, *space_scalar));
	DataSet	*dataset_vol	= new DataSet(file.createDataSet("Volume", PredType::NATIVE_FLOAT, *space_scalar));
	DataSet	*dataset_pos	= new DataSet(file.createDataSet("Position", PredType::NATIVE_FLOAT, *space_vector));
	DataSet	*dataset_vel	= new DataSet(file.createDataSet("Velocity", PredType::NATIVE_FLOAT, *space_vector));
	DataSet	*dataset_psize	= new DataSet(file.createDataSet("Psize", PredType::NATIVE_FLOAT, *space_vector));
	DataSet	*dataset_s		= new DataSet(file.createDataSet("Stress", PredType::NATIVE_FLOAT, *space_tensor));
	DataSet	*dataset_td		= new DataSet(file.createDataSet("Derformation_Tensor", PredType::NATIVE_FLOAT, *space_true_tensor));

	DataSet	*dataset_pre	= new DataSet(file.createDataSet("Pressure", PredType::NATIVE_FLOAT, *space_scalar));

	dataset_tag->write(tag_h5, PredType::NATIVE_FLOAT);
	dataset_m->write(m_h5, PredType::NATIVE_FLOAT);
	dataset_vol->write(vol_h5, PredType::NATIVE_FLOAT);
	dataset_pos->write(pos_h5, PredType::NATIVE_FLOAT);
	dataset_vel->write(vel_h5, PredType::NATIVE_FLOAT);
	dataset_psize->write(psize_h5, PredType::NATIVE_FLOAT);
	dataset_s->write(s_h5, PredType::NATIVE_FLOAT);
	dataset_td->write(td_h5, PredType::NATIVE_FLOAT);

	dataset_pre->write(pressure_h5, PredType::NATIVE_FLOAT);

	delete space_scalar;
	delete space_vector;
	delete space_tensor;
	delete space_true_tensor;
	delete dataset_tag;
	delete dataset_m;
	delete dataset_vol;
	delete dataset_pos;
	delete dataset_vel;
	delete dataset_psize;
	delete dataset_s;
	delete dataset_td;

	delete dataset_pre;

	delete[] tag_h5;
	delete[] m_h5;
	delete[] vol_h5;
	delete[] pos_h5;
	delete[] vel_h5;
	delete[] psize_h5;
	delete[] s_h5;
	delete[] td_h5;

	delete[] pressure_h5;

	file.close();

	std::string file_name_xmf = name+"_MPM_PARTICLE_"+out.str()+".xmf";

    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Lp.size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Lp.size() << " 3\" >\n";
    oss << "        " << file_name_h5 <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Vol\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Vol \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lp.size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Stress\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}

template<int SType, int D>
inline void MPM<SType, D>::WriteFileH5Node(string name, int n)
{
    std::stringstream   out;                    //convert int to string for file name.
    out << std::setw(6) << std::setfill('0') << n;          
    std::string file_name_h5 = name+"_MPM_NODE_"+out.str()+".h5";

	H5File  file(file_name_h5, H5F_ACC_TRUNC);
    hsize_t nx = Nx+1;
    hsize_t ny = Ny+1;
    hsize_t nz = Nz+1;

    float dx = Dx;
    float ox = Origin(0);
    float oy = Origin(1);
    float oz = Origin(2);

    hsize_t dims_scalar[3] = {nx, ny, nz};          //create data space.
    hsize_t dims_vector[4] = {nx, ny, nz, 3};       //create data space.

    int rank_scalar = sizeof(dims_scalar) / sizeof(hsize_t);
    int rank_vector = sizeof(dims_vector) / sizeof(hsize_t);

    DataSpace   *space_scalar = new DataSpace(rank_scalar, dims_scalar);
    DataSpace   *space_vector = new DataSpace(rank_vector, dims_vector);

    DataSet *dataset_m = new DataSet(file.createDataSet("Mass", PredType::NATIVE_FLOAT, *space_scalar));
	DataSet *dataset_sdf = new DataSet(file.createDataSet("SDF", PredType::NATIVE_FLOAT, *space_scalar));
    DataSet *dataset_vel = new DataSet(file.createDataSet("Velocity", PredType::NATIVE_FLOAT, *space_vector));

    float* mass_h5  = new float[Nnode];
	float* sdf_h5   = new float[Nnode];
    float* vel_h5   = new float[3*Nnode];

	for (size_t l=0; l<Nnode; ++l) {		
		mass_h5[l]      = Ln[l]->M;
		sdf_h5[l]		= Ln[l]->SDF;
		vel_h5[3*l]     = Ln[l]->V[0];
		vel_h5[3*l+1]   = Ln[l]->V[1];
		vel_h5[3*l+2]   = Ln[l]->V[2];
	}
    dataset_m->write(mass_h5, PredType::NATIVE_FLOAT);
	dataset_sdf->write(sdf_h5, PredType::NATIVE_FLOAT);
    dataset_vel->write(vel_h5, PredType::NATIVE_FLOAT);

    delete space_scalar;
    delete space_vector;
    delete dataset_m;
	delete dataset_sdf;
    delete dataset_vel;

    delete[] mass_h5;
	delete[] sdf_h5;
    delete[] vel_h5;

    file.close();

    std::string file_name_xmf = name+"_MPM_NODE_"+out.str()+".xmf";
    std::ofstream oss;
    oss.open(file_name_xmf);
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM_NODE\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << nz << " " << ny << " " << nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << ox << " " << oy  << " " << oz  << "\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << dx << " " << dx  << " " << dx  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Mass\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Mass \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"SDF\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/SDF \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << nz << " " << ny << " " << nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << file_name_h5 <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    oss.close();
}


// inline void MPM::LoadMPMFromH5(string fname, double ratio)
// {
// 	cout << "========= Start loading MPM particles from " << fname << "==============" << endl;
// 	H5std_string FILE_NAME( fname );
// 	H5std_string DATASET_NAME_POS( "Position" );
// 	H5File file_pos( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_pos = file_pos.openDataSet( DATASET_NAME_POS );
// 	DataSpace dataspace_pos = dataset_pos.getSpace();
//     hsize_t dims_pos[2];
//     dataspace_pos.getSimpleExtentDims( dims_pos, NULL);
//     hsize_t dimsm_pos = dims_pos[0];
//     cout <<"Position" << endl;

// 	H5std_string DATASET_NAME_VEL( "Velocity" );
// 	H5File file_vel( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_vel = file_pos.openDataSet( DATASET_NAME_VEL );
// 	DataSpace dataspace_vel = dataset_vel.getSpace();
//     hsize_t dims_vel[2];
//     dataspace_vel.getSimpleExtentDims( dims_vel, NULL);
//     hsize_t dimsm_vel = dims_vel[0];
//     cout <<"Velocity" << endl;

// 	H5std_string DATASET_NAME_PSIZE( "Psize" );
// 	H5File file_psize( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_psize = file_pos.openDataSet( DATASET_NAME_PSIZE );
// 	DataSpace dataspace_psize = dataset_psize.getSpace();
//     hsize_t dims_psize[2];
//     dataspace_psize.getSimpleExtentDims( dims_psize, NULL);
//     hsize_t dimsm_psize = dims_psize[0];
//     cout <<"Psize" << endl;

// 	H5std_string DATASET_NAME_S( "Stress" );
// 	H5File file_s( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_s = file_pos.openDataSet( DATASET_NAME_S );
// 	DataSpace dataspace_s = dataset_s.getSpace();
//     hsize_t dims_s[2];
//     dataspace_s.getSimpleExtentDims( dims_s, NULL);
//     hsize_t dimsm_s = dims_s[0];
//     cout <<"Stress" << endl;

// 	H5std_string DATASET_NAME_TD( "Derformation_Tensor" );
// 	H5File file_td( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_td = file_pos.openDataSet( DATASET_NAME_TD );
// 	DataSpace dataspace_td = dataset_td.getSpace();
//     hsize_t dims_td[2];
//     dataspace_td.getSimpleExtentDims( dims_td, NULL);
//     hsize_t dimsm_td = dims_td[0];
//     cout <<"Derformation_Tensor" << endl;

// 	H5std_string DATASET_NAME_M( "Mass" );
// 	H5File file_m( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_m = file_m.openDataSet( DATASET_NAME_M );
// 	DataSpace dataspace_m = dataset_m.getSpace();
//     hsize_t dims_m[2];
//     dataspace_m.getSimpleExtentDims( dims_m, NULL);
//     hsize_t dimsm_m = dims_m[0];
//     cout <<"Mass" << endl;

// 	H5std_string DATASET_NAME_VOL( "Volume" );
// 	H5File file_vol( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_vol = file_vol.openDataSet( DATASET_NAME_VOL );
// 	DataSpace dataspace_vol = dataset_vol.getSpace();
//     hsize_t dims_vol[2];
//     dataspace_vol.getSimpleExtentDims( dims_vol, NULL);
//     hsize_t dimsm_vol = dims_vol[0];
//     cout <<"Volume" << endl;

// 	H5std_string DATASET_NAME_TAG( "Tag" );
// 	H5File file_tag( FILE_NAME, H5F_ACC_RDONLY );
// 	DataSet dataset_tag = file_tag.openDataSet( DATASET_NAME_TAG );
// 	DataSpace dataspace_tag = dataset_tag.getSpace();
//     hsize_t dims_tag[2];
//     dataspace_tag.getSimpleExtentDims( dims_tag, NULL);
//     hsize_t dimsm_tag = dims_tag[0];
//     cout <<"Tag" << endl;

//     double* data_pos = new double[dimsm_pos];
//     dataset_pos.read( data_pos, PredType::NATIVE_FLOAT, dataspace_pos, dataspace_pos );
//     cout <<"data_pos" << endl;

//     double* data_vel = new double[dimsm_vel];
//     dataset_vel.read( data_vel, PredType::NATIVE_FLOAT, dataspace_vel, dataspace_vel );
//     cout <<"data_vel" << endl;

//     double* data_s = new double[dimsm_s];
//     dataset_s.read( data_s, PredType::NATIVE_FLOAT, dataspace_s, dataspace_s );
//     cout <<"data_s" << endl;

//     double* data_td = new double[dimsm_td];
//     dataset_td.read( data_td, PredType::NATIVE_FLOAT, dataspace_td, dataspace_td );
//     cout <<"data_td" << endl;

//     double* data_m = new double[dimsm_m];
//     dataset_m.read( data_m, PredType::NATIVE_FLOAT, dataspace_m, dataspace_m );
//     cout <<"data_m" << endl;

//     double* data_tag = new double[dimsm_tag];
//     dataset_tag.read( data_tag, PredType::NATIVE_FLOAT, dataspace_tag, dataspace_tag );
//     cout <<"data_tag" << endl;

//     double* data_vol = new double[dimsm_vol];
//     dataset_vol.read( data_vol, PredType::NATIVE_FLOAT, dataspace_vol, dataspace_vol );
//     cout <<"data_vol" << endl;

//     double* data_psize = new double[dimsm_psize];
//     dataset_psize.read( data_psize, PredType::NATIVE_FLOAT, dataspace_psize, dataspace_psize );
//     cout <<"data_psize" << endl;

//     int np = dimsm_pos/3;
//     for (int i=0; i<np; ++i)
//     {
//     	Vector3d pos (data_pos[3*i], data_pos[3*i+1], data_pos[3*i+2]);
//     	Vector3d vel (data_vel[3*i], data_vel[3*i+1], data_vel[3*i+2]);
//     	Vector3d psize (data_psize[3*i], data_psize[3*i+1], data_psize[3*i+2]);
//     	Matrix3d stress;
// 		stress(0,0) = data_s[6*i];
// 		stress(0,1) = stress(1,0) = data_s[6*i+1];
// 		stress(0,2) = stress(2,0) = data_s[6*i+2];
// 		stress(1,1) = data_s[6*i+3];
// 		stress(2,1) = stress(1,2) = data_s[6*i+4];
// 		stress(2,2) = data_s[6*i+5];
//     	Matrix3d td;
// 		td(0,0) = data_td[9*i];
// 		td(0,1) = data_td[9*i+1];
// 		td(0,2) = data_td[9*i+2];
// 		td(1,0) = data_td[9*i+3];
// 		td(1,1) = data_td[9*i+4];
// 		td(1,2) = data_td[9*i+5];
// 		td(2,0) = data_td[9*i+6];
// 		td(2,1) = data_td[9*i+7];
// 		td(2,2) = data_td[9*i+8];

//     	double m = data_m[i];
//     	double vol = data_vol[i];
//     	int tag = (int) data_tag[i];
//     	AddParticle(tag, pos, m);
//     	Lp[Lp.size()-1]->V = vel;
//     	Lp[Lp.size()-1]->PSize = psize;
//     	Lp[Lp.size()-1]->Vol = vol;
//     	Lp[Lp.size()-1]->Stress = stress;
//     	Lp[Lp.size()-1]->Td = td;
//     	Lp[Lp.size()-1]->Vol0 = 1.;
//     	for (size_t d=0; d<D; ++d)
//     	{
//     		// Lp[Lp.size()-1]->Vol0 *= ratio;
//     		Lp[Lp.size()-1]->PSize0(d) = 0.5*ratio;
//     	}
//     }

//     delete data_pos;
//     delete data_vel;
//     delete data_s;
//     delete data_td;
//     delete data_m;
//     delete data_tag;
//     delete data_vol;
//     delete data_psize;

//     cout << "========= Loaded "<< Lp.size()<< " MPM particles from " << fname << "==============" << endl;
// }