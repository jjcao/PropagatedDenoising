#include"ShortestPropagationMeshFiltering.h"
#include<iostream>
//#include <limits>
//#include <algorithm>
//#include "Eigen/Dense"
//#include "Eigen/Sparse"

void ShortestPropagationMeshFiltering::updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals)
{
	filtered_normals.resize((int)mesh.n_faces());
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


	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	FaceNeighborType face_neighbor_type = static_cast<FaceNeighborType>(face_neighbor_index);
	double radius;
	if (face_neighbor_type == kRadiusBased)
		radius = getAveragefaceCenterDistance(mesh) * multiple_radius;

	std::vector<std::vector<TriMesh::FaceHandle> > all_face_neighbor((int)mesh.n_faces());
	getAllFaceNeighbor(mesh, all_face_neighbor, face_neighbor_type, include_central_face, radius);
	getFaceNormal(mesh, filtered_normals);

	std::vector<double> face_area((int)mesh.n_faces());
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	std::vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());//上一次的法线

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	//Eigen::SparseMatrix<std::vector<TriMesh::FaceHandle>> path_matrix((int)mesh.n_faces(), (int)mesh.n_faces());
	//std::vector<TriMesh::FaceHandle> tmpPath = path_matrix.coeff(0, 0);
	//getAllFaceNeighborAdjacencyList(mesh, all_face_neighbor, all_face_neighbor_adja_list);
	for (int iter = 0; iter < normal_iteration_number; ++iter)
	{
		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		getFaceNormal(mesh, previous_normals);

		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			int index = f_it->idx();
			std::vector<TriMesh::FaceHandle> face_neighbor = all_face_neighbor[index];

			//---------------------------------------得到局部图-------------------------------------
			AdjacencyList localGraph;
			std::vector<int> local2GlobalIdx;
			std::map<int, int> global2localIdx;
			buildLocalGraph(mesh, face_neighbor, localGraph, local2GlobalIdx, global2localIdx);
			Dijkstra dij;
			std::vector<double> distance;
			dij.computeDistances(localGraph, global2localIdx[index], distance);
			//--------------------------------------------------------------------------------------

			TriMesh::Normal filtered_normal(0.0, 0.0, 0.0);
			for (int j = 0; j < (int)face_neighbor.size(); ++j)
			{
				int currentIdx = face_neighbor[j].idx();
				//TriMesh::Normal direct_center = (face_centroid[index] - face_centroid[currentIdx]).normalize();
				std::vector<int> facePath;
				double weight = 0.0;
				if (face_neighbor[j] != *f_it)
				{	
					dij.computePath(localGraph, global2localIdx[currentIdx], facePath);

					double sumPF1 = 0.0, sumPF2 = 0.0;
					//前后两项的法线差异
					for (int k = (int)facePath.size() - 1; k > 0; --k)
					{
						int preIndex = facePath[k];
						int nexIndex = facePath[k - 1];
						double temp11 = (previous_normals[nexIndex] - previous_normals[preIndex]).length();
						sumPF1 += temp11*temp11;
					}
					//累积路径差异
					for (int k = (int)facePath.size() - 1; k >= 0; k--)
					{
						int currentIndex = facePath[k];
						double temp22 = (previous_normals[currentIndex] - previous_normals[index]).length();
						sumPF2 += temp22*temp22;
					}

					//double sigmaA = 1, sigmaR = 1;
					weight = GaussianWeight(sqrt(sumPF1), sigma_s) * GaussianWeight(sqrt(sumPF2), sigma_r);
				}
				else
				{
					weight = 1.0;
				}
				filtered_normal += weight * face_area[currentIdx] * previous_normals[currentIdx];
			}
			if (face_neighbor.size())
				filtered_normals[index] = filtered_normal.normalize();
		}
		// immediate update vertex position
		updateVertexPosition(mesh, filtered_normals, vertex_iteration_number, false);
	}
}

void ShortestPropagationMeshFiltering::updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals)
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
void ShortestPropagationMeshFiltering::denoise()
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
ShortestPropagationMeshFiltering::ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set)
	: MeshDenoisingBase(_data_manager, _parameter_set)
{
	initParameters();
}
void ShortestPropagationMeshFiltering::initParameters()
{
	parameter_set_->removeAllParameter();
	QStringList strList_DenoiseType;
	strList_DenoiseType.push_back(QString("Local"));
	strList_DenoiseType.push_back(QString("Global"));

	parameter_set_->addParameter(QString("Denoise Type"), strList_DenoiseType, 0, QString("Denoise Type"), QString("The type of denoise method."));

	QStringList strList_FaceNeighborType;
	strList_FaceNeighborType.push_back(QString("geometrical"));
	strList_FaceNeighborType.push_back(QString("topological"));

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

double ShortestPropagationMeshFiltering::GaussianWeight(double distance, double sigma)
{
	return std::exp(-0.5 * distance * distance / (sigma * sigma));
}

double ShortestPropagationMeshFiltering::NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2)
{
	return (n1 - n2).length();
}

//网格上所有面之间的面心距离的平均值的倍数，作为衡量基准，
//一般取一倍，所以值就是面心距离的平均值
double ShortestPropagationMeshFiltering::getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh)
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

void ShortestPropagationMeshFiltering::buildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, 
	AdjacencyList &localGraph, std::vector<int>& local2GlobalIdx, std::map<int, int>& global2localIdx)
{
	local2GlobalIdx.clear(); global2localIdx.clear();
	localGraph.clear(); 
	localGraph.resize(face_neighbor.size());

	int i(0);
	for (std::vector<TriMesh::FaceHandle>::iterator iter = face_neighbor.begin(); iter != face_neighbor.end(); ++iter, ++i)
	{
		local2GlobalIdx.push_back(iter->idx());
		global2localIdx.insert(std::make_pair(iter->idx(),i));
	}
		
	i=0;
	for (std::vector<TriMesh::FaceHandle>::iterator iter = face_neighbor.begin(); iter != face_neighbor.end(); ++iter, ++i)
	{
		TriMesh::Point c1 = mesh.calc_face_centroid(*iter);
		std::vector<TriMesh::FaceHandle> locFaceNeighbor;
		// todo edgeBased v.s. vertexBased
		getFaceNeighbor(mesh, *iter, kEdgeBased, locFaceNeighbor);
		
		for (std::vector<TriMesh::FaceHandle>::iterator it = locFaceNeighbor.begin(); it != locFaceNeighbor.end(); ++it)
		{
			TriMesh::Point c2 = mesh.calc_face_centroid(*it);
			localGraph[i].push_back(std::make_pair(global2localIdx[it->idx()], (c1 - c2).length())); 
		}
	}
}

