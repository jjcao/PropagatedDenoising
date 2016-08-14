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
	double radius = multiple_radius;
	FaceNeighborType face_neighbor_type = static_cast<FaceNeighborType>(face_neighbor_index);
	if (face_neighbor_type == kRadiusBased)
		radius = getAveragefaceCenterDistance(mesh) * multiple_radius;

	setAllFaceNeighbor(mesh, face_neighbor_type, include_central_face, radius);
	
	getFaceNormal(mesh, filtered_normals);

	std::vector<double> face_area((int)mesh.n_faces());
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	std::vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());//上一次的法线
	std::vector<TriMesh::Normal> guided_normals((int)mesh.n_faces());
	std::vector<std::pair<double, TriMesh::Normal> > range_and_mean_normal((int)mesh.n_faces());

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////
	for (int iter = 0; iter < normal_iteration_number; ++iter)
	{
		std::cout << "iteration: " << iter << std::endl;

		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		getFaceNormal(mesh, previous_normals);

		getGuidedNormals(mesh, face_area, previous_normals, range_and_mean_normal, guided_normals);

		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			int index = f_it->idx();		
			
			std::vector<std::vector<int> > facePaths;
			computeGlobalPath(mesh, f_it, face_centroid, guided_normals, facePaths);
			//sigma_r = calculateSigma(previous_normals, f_it, iter, sigma_s);


			TriMesh::Normal filteredNormal(0.0, 0.0, 0.0);
			//std::vector<double> sumAs;
			//std::vector<double> sumRs;
			//std::vector<int> neighbors;
			for (int j = 0; j < facePaths.size(); ++j)
			{
				std::vector<int> &facePath = facePaths[j];
				int neighborIdx = facePath[0];//remember facePath[0] is neighbor!
				
				double weight = 0.0;
				double sumA = 0.0;// adjacent photometric relationship
				double sumR = 0.0;// photometric relationship
				for (int k = (int)facePath.size() - 1; k > 0; --k)
				{
					int preIndex = facePath[k];
					int nexIndex = facePath[k - 1];
					double distD = (guided_normals[nexIndex] - guided_normals[preIndex]).length();
					sumA += distD*distD;
					double distR = (guided_normals[preIndex] - guided_normals[index]).length();
					sumR += distR*distR;
				}

				//sumAs.push_back(sumA);
				//sumRs.push_back(sumR);
				//neighbors.push_back(neighborIdx);

				weight = GaussianWeight(sqrt(sumA), sigma_r) * GaussianWeight(sqrt(sumR), sigma_r);

				//--------why use previous_normals?
				filteredNormal += weight * face_area[neighborIdx] * previous_normals[neighborIdx];
			}
			//double sigma1 = calculateSigmaSandR(sumAs,sigma_s);
			//double sigma2 = calculateSigmaSandR(sumRs,sigma_s);
			//for (int neighbor_num = 0; neighbor_num < neighbors.size(); ++neighbor_num)
			//{
			//	double weight = GaussianWeight(sqrt(sumAs[neighbor_num]), sigma1) * GaussianWeight(sqrt(sumRs[neighbor_num]), sigma2);
			//	filteredNormal += weight * face_area[neighbors[neighbor_num]] * previous_normals[neighbors[neighbor_num]];
			//}

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

	parameter_set_->addParameter(QString("Multiple(* avg face dis.)"), 3.0, QString("Multiple(* avg face dis.)"), QString("Radius for search geometrical neighbor of the face."),
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
double PropagationMeshFiltering::calculateSigmaSandR(std::vector<double> &datas, double smoothness)
{
	double aver = 0.0;
	int len = datas.size();
	for (int i = 0; i < len; ++i)
	{
		aver += datas[i];
	}
	aver = aver / len;

	double sigma = 0.0;
	for (int i = 0; i < len; ++i)
	{
		sigma += (datas[i] - aver) * (datas[i] - aver);
	}

	return sqrt(sigma / len) + smoothness;
}

void PropagationMeshFiltering::getFaceNeighborInnerEdge(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &face_neighbor, std::vector<TriMesh::EdgeHandle> &inner_edge)
{
	inner_edge.clear();
	std::vector<bool> edge_flag((int)mesh.n_edges(), false);
	std::vector<bool> face_flag((int)mesh.n_faces(), false);

	for (int i = 0; i < (int)face_neighbor.size(); i++)
		face_flag[face_neighbor[i].idx()] = true;

	for (int i = 0; i < (int)face_neighbor.size(); i++)
	{
		for (TriMesh::FaceEdgeIter fe_it = mesh.fe_iter(face_neighbor[i]); fe_it.is_valid(); fe_it++)
		{
			if ((!edge_flag[fe_it->idx()]) && (!mesh.is_boundary(*fe_it)))
			{
				edge_flag[fe_it->idx()] = true;
				TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(*fe_it, 0);
				TriMesh::FaceHandle f = mesh.face_handle(heh);
				TriMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
				TriMesh::FaceHandle fo = mesh.face_handle(heho);
				if (face_flag[f.idx()] && face_flag[fo.idx()])
					inner_edge.push_back(*fe_it);
			}
		}
	}
}

void PropagationMeshFiltering::getRangeAndMeanNormal(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_guided_neighbor,
	std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
	std::vector<std::pair<double, TriMesh::Normal> > &range_and_mean_normal)
{
	const double epsilon = 1.0e-9;

	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		int index = f_it->idx();
		std::vector<TriMesh::FaceHandle> face_neighbor = all_guided_neighbor[index];

		double metric = 0.0;
		TriMesh::Normal average_normal(0.0, 0.0, 0.0);
		double maxdiff = -1.0;

		for (int i = 0; i < (int)face_neighbor.size(); i++)
		{
			int index_i = face_neighbor[i].idx();
			double area_weight = face_area[index_i];
			TriMesh::Normal ni = normals[index_i];
			average_normal += ni * area_weight;

			for (int j = i + 1; j < (int)face_neighbor.size(); j++)
			{
				int index_j = face_neighbor[j].idx();
				TriMesh::Normal nj = normals[index_j];
				double diff = NormalDistance(ni, nj);

				if (diff > maxdiff)
				{
					maxdiff = diff;
				}
			}
		}

		std::vector<TriMesh::EdgeHandle> inner_edge_handle;
		getFaceNeighborInnerEdge(mesh, face_neighbor, inner_edge_handle);
		double sum_tv = 0.0, max_tv = -1.0;
		for (int i = 0; i < (int)inner_edge_handle.size(); i++)
		{
			TriMesh::HalfedgeHandle heh = mesh.halfedge_handle(inner_edge_handle[i], 0);
			TriMesh::FaceHandle f = mesh.face_handle(heh);
			TriMesh::Normal n1 = normals[f.idx()];
			TriMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
			TriMesh::FaceHandle fo = mesh.face_handle(heho);
			TriMesh::Normal n2 = normals[fo.idx()];
			double current_tv = NormalDistance(n1, n2);
			max_tv = (current_tv > max_tv) ? current_tv : max_tv;
			sum_tv += current_tv;
		}

		average_normal.normalize();
		metric = maxdiff * max_tv / (sum_tv + epsilon);

		range_and_mean_normal[index] = std::make_pair(metric, average_normal);
	}
}

