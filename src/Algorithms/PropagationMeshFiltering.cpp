#include "PropagationMeshFiltering.h"
#include<iostream>
#include <sstream>

void PropagationMeshFiltering::updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals)
{
	// get parameter for local scheme normal update
	int face_neighbor_index;
	if (!parameter_set_->getStringListIndex(QString("Face Neighbor"), face_neighbor_index))
		return;
	bool include_central_face;
	if (!parameter_set_->getValue(QString("include central face"), include_central_face))
		return;
	double multiple_radius;
	if (!parameter_set_->getValue(QString("Multiple(* avg face dis.)"), multiple_radius))
		return;
	double sigma_s;
	if (!parameter_set_->getValue(QString("Multiple(* sigma_s)"), sigma_s)) // use it as sigma_s actually.
		return;
	int normal_iteration_number;
	if (!parameter_set_->getValue(QString("(Local)Normal Iteration Num."), normal_iteration_number))
		return;
	double sigma_r;
	if (!parameter_set_->getValue(QString("sigma_r"), sigma_r))
		return;
	int vertex_iteration_number;
	if (!parameter_set_->getValue(QString("Vertex Iteration Num."), vertex_iteration_number))
		return;

	std::cout << "PropagationMeshFiltering with Multiple(* avg face dis.) = " << multiple_radius << std::endl;

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	FaceNeighborType face_neighbor_type = static_cast<FaceNeighborType>(face_neighbor_index);
	double radius = multiple_radius;

	//FaceNeighborType face_neighbor_type = face_neighbor_index == 0 ? kRadiusBased : kVertexBased;

	//double radius;
	//if (face_neighbor_type == kRadiusBased)
	//	radius = getRadius(multiple_radius, mesh);

	if (face_neighbor_type == kRadiusBased)
		radius = getAveragefaceCenterDistance(mesh) * multiple_radius;
	setAllFaceNeighbor(mesh, face_neighbor_type, include_central_face, radius);
	
	getFaceNormal(mesh, filtered_normals);
	std::vector<double> face_area((int)mesh.n_faces());
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	std::vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());//上一次的法线

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	for (int iter = 0; iter < normal_iteration_number; ++iter)
	{
		std::cout << "iteration: " << iter << std::endl;

		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		getFaceNormal(mesh, previous_normals);

		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			int index = f_it->idx();		
			
			std::vector<std::vector<int> > facePaths;
			computeGlobalPath(mesh, f_it, face_centroid, previous_normals, facePaths);
			sigma_r = calculateSigma(previous_normals, f_it, iter, sigma_s);


			TriMesh::Normal filteredNormal(0.0, 0.0, 0.0);
			for (int j = 0; j < facePaths.size(); ++j)
			{
				std::vector<int> &facePath = facePaths[j];
				int neighborIdx = facePath[0];//remember facePath[0] is neighbor!

				double weight = 0.0;
				if (index != neighborIdx)
				{
					double sumA = 0.0;// adjacent photometric relationship
					double sumR = 0.0;// photometric relationship
					for (int k = (int)facePath.size() - 1; k > 0; --k)
					{
						int preIndex = facePath[k];
						int nexIndex = facePath[k - 1];
						double distD = (previous_normals[nexIndex] - previous_normals[preIndex]).length();
						sumA += distD*distD;
						double distR = (previous_normals[preIndex] - previous_normals[index]).length();
						sumR += distR*distR;
					}

					//得到SigmaS 和SigmaR(//SigmaR指的是法向差的范围，这个没有给出某种统计值，指导的这篇默认用的是0.35	//)
					// todo: sigma for sumA should be larger than sumR?
					weight = GaussianWeight(sqrt(sumA), sigma_r) * GaussianWeight(sqrt(sumR), sigma_s*sigma_r);
				}
				else{
					neighborIdx = index;
				}

				filteredNormal += weight * face_area[neighborIdx] * previous_normals[neighborIdx];
			}

			filtered_normals[index] = filteredNormal.normalize();
		}

		// immediate update vertex position
		updateVertexPosition(mesh, filtered_normals, vertex_iteration_number, false);

		//checkBadFace(mesh);//对坏的面进行微小的扰动

		// for debug
		if (iter % 4 == 0)
		{
			data_manager_->setMesh(mesh);
			std::stringstream ss;  ss << iter << ".obj";
			data_manager_->ExportMeshToFile(ss.str());
		}
	}
}

void PropagationMeshFiltering::updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals)
{
	// get parameter for normal update
	int denoise_index;
	if (!parameter_set_->getStringListIndex(QString("Denoise Type"), denoise_index))
		return;

	DenoiseType denoise_type = (denoise_index == 0) ? kLocal : kGlobal;

	if (denoise_type == kLocal)
		updateFilteredNormalsLocalScheme(mesh, filtered_normals);
	//else if (denoise_type == kGlobal)
	//	updateFilteredNormalsGlobalScheme(mesh, filtered_normals);

}
void PropagationMeshFiltering::denoise()
{
	// get data
	TriMesh mesh = data_manager_->getNoisyMesh();

	if (mesh.n_vertices() == 0)
		return;

	// update face normal
	std::vector<TriMesh::Normal> filtered_normals;
	updateFilteredNormals(mesh, filtered_normals);

	// update data
	data_manager_->setMesh(mesh);
	data_manager_->setDenoisedMesh(mesh);
}


PropagationMeshFiltering::PropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set) : MeshDenoisingBase(_data_manager, _parameter_set)
{
	initParameters();
}

void PropagationMeshFiltering::initParameters()
{
	parameter_set_->removeAllParameter();
	QStringList strList_DenoiseType;
	strList_DenoiseType.push_back(QString("Local"));
	strList_DenoiseType.push_back(QString("Global"));

	parameter_set_->addParameter(QString("Denoise Type"), strList_DenoiseType, 0, QString("Denoise Type"), QString("The type of denoise method."));

	QStringList strList_FaceNeighborType;
	strList_FaceNeighborType.push_back(QString("kRadiusBased"));
	strList_FaceNeighborType.push_back(QString("kVertexBased"));
	strList_FaceNeighborType.push_back(QString("kEdgeBased"));
	strList_FaceNeighborType.push_back(QString("kFaceRingBased"));

	parameter_set_->addParameter(QString("Face Neighbor"), strList_FaceNeighborType, 0, QString("Face Neighbor"), QString("The type of the neighbor of the face."));
	parameter_set_->addParameter(QString("include central face"), true, QString("include central face"), QString("Include the central face of the neighbor or not."));

	parameter_set_->addParameter(QString("Multiple(* avg face dis.)"), 2.0, QString("Multiple(* avg face dis.)"), QString("Radius for search geometrical neighbor of the face."),
		true, 0.1, 10.0);
	parameter_set_->addParameter(QString("Multiple(* sigma_s)"), 1.0, QString("Multiple(* sigma_s)"), QString("Standard deviation of spatial weight."),
		true, 1.0e-9, 10.0);

	parameter_set_->addParameter(QString("sigma_r"), 0.35, QString("sigma_r"), QString("Standard deviation of range weight."),
		true, 1.0e-9, 10.0);

	parameter_set_->addParameter(QString("(Local)Normal Iteration Num."), 20, QString("(Local)Normal Iteration Num."), QString("The iteration number of local scheme for filtering face normal."),
		true, 1, 500);
	parameter_set_->addParameter(QString("(Global)Normal Iteration Num."), 1, QString("(Global)Normal Iteration Num."), QString("The iteration number of global scheme for filtering face normal."),
		true, 1, 50);

	parameter_set_->addParameter(QString("smoothness"), 0.01, QString("smoothness"), QString("Smoothness parameter for global filtering."),
		true, 1.0e-9, 1.0);

	parameter_set_->addParameter(QString("Vertex Iteration Num."), 10, QString("Vertex Iteration Num."), QString("The iteration number for filtering vertex position."),
		true, 1, 500);


	parameter_set_->setName(QString("Propagation Mesh Filtering"));
	parameter_set_->setLabel(QString("Propagation Mesh Filtering"));
	parameter_set_->setIntroduction(QString("Propagation Mesh Filtering -- Parameters"));
}

double PropagationMeshFiltering::GaussianWeight(double distance, double sigma)
{
	return std::exp(-0.5 * distance * distance / (sigma * sigma));
}

double PropagationMeshFiltering::NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2)
{
	return (n1 - n2).length();
}

double PropagationMeshFiltering::getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh)
{

	double sigma_s = 0.0, num = 0.0;
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point fi = centroid[f_it->idx()];
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
		{
			TriMesh::Point fj = centroid[ff_it->idx()];
			sigma_s += (fj - fi).length();
			num++;
		}
	}
	return sigma_s * multiple / num;
}