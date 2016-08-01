#include"ShortestPropagationMeshFiltering.h"


void ShortestPropagationMeshFiltering::initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, int sourceIdxGlobal)
{
	_local2GlobalIdx.clear(); _global2localIdx.clear();
	_localGraph.clear();
	_localGraph.resize(faceNeighbor.size());


	int i(0);
	for (std::vector<TriMesh::FaceHandle>::iterator iter = faceNeighbor.begin(); iter != faceNeighbor.end(); ++iter, ++i)
	{
		_local2GlobalIdx.push_back(iter->idx());
		_global2localIdx.insert(std::make_pair(iter->idx(), i));
	}

	i = 0;
	for (std::vector<TriMesh::FaceHandle>::iterator iter = faceNeighbor.begin(); iter != faceNeighbor.end(); ++iter, ++i)
	{
		TriMesh::Point c1 = mesh.calc_face_centroid(*iter);
		std::vector<TriMesh::FaceHandle> locFaceNeighbor;
		// todo edgeBased v.s. vertexBased
		getFaceNeighbor(mesh, *iter, kEdgeBased, locFaceNeighbor);

		for (std::vector<TriMesh::FaceHandle>::iterator it = locFaceNeighbor.begin(); it != locFaceNeighbor.end(); ++it)
		{
			TriMesh::Point c2 = mesh.calc_face_centroid(*it);
			_localGraph[i].push_back(std::make_pair(_global2localIdx[it->idx()], (c1 - c2).length()));
		}
	}

	std::vector<double> distance;
	_dij.computeDistances(_localGraph, _global2localIdx[sourceIdxGlobal], distance);
}

void ShortestPropagationMeshFiltering::computeGlobalPath(int targetIdxGlobal, std::vector<int>& facePath)
{
	_dij.computePath(_localGraph, _global2localIdx[targetIdxGlobal], facePath);
	for (int k = 0; k < facePath.size(); ++k)
		facePath[k] = _local2GlobalIdx[facePath[k]];
}

double ShortestPropagationMeshFiltering::calculateSigma(std::vector<TriMesh::Normal> &face_normals, std::vector<TriMesh::FaceHandle> &faceNeighbor, double smoothness)
{
	//自适应局部区域的法线情况，计算法线的方差，用方差来作为高斯的方差，（注高斯函数中的sigma是标准差）
	TriMesh::Normal aver_local_normal(0.0, 0.0, 0.0);
	int len = faceNeighbor.size();
	for (int st = 0; st < len; st++)
	{
		aver_local_normal += face_normals[faceNeighbor[st].idx()];
	}
	aver_local_normal = aver_local_normal / len;

	double stdard = 0.0;
	for (int st = 0; st < len; st++)
	{
		double dtemp = (aver_local_normal - face_normals[faceNeighbor[st].idx()]).length();
		stdard += dtemp * dtemp;
		//stdard += NormalDistance(aver_local_normal, previous_normals[projections[i][st].face_index]);
	}
	return sqrt(stdard / len) + smoothness;//sigma_s替代， 作为光滑的参数
}

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

	std::cout << "PropagationMeshFiltering with Multiple(* avg face dis.) = " << multiple_radius << std::endl;

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
		std::cout << "iteration: " << iter << std::endl;

		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		getFaceNormal(mesh, previous_normals);

		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			int index = f_it->idx();
			std::vector<TriMesh::FaceHandle> faceNeighbor = all_face_neighbor[index];

			initLocalGraph(mesh, faceNeighbor, index);
			double sigma_r = calculateSigma(previous_normals, faceNeighbor, sigma_s);

			TriMesh::Normal filteredNormal(0.0, 0.0, 0.0);
			for (int j = 0; j < (int)faceNeighbor.size(); ++j)
			{
				int currentIdx = faceNeighbor[j].idx();
				//TriMesh::Normal direct_center = (face_centroid[index] - face_centroid[currentIdx]).normalize();
				std::vector<int> facePath;
				double weight = 0.0;
				if (faceNeighbor[j] != *f_it)
				{
					computeGlobalPath(currentIdx, facePath);

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
					for (int k = (int)facePath.size() - 1; k >= 0; --k)
					{
						int currentIndex = facePath[k];
						double temp22 = (previous_normals[currentIndex] - previous_normals[index]).length();
						sumPF2 += temp22*temp22;
					}

					//double sigmaA = 1, sigmaR = 1;
					//得到SigmaS 和SigmaR(//SigmaR指的是法向差的范围，这个没有给出某种统计值，指导的这篇默认用的是0.35	//)
					weight = GaussianWeight(sqrt(sumPF1), sigma_r) * GaussianWeight(sqrt(sumPF2), sigma_r);
				}
				else
				{
					weight = 1.0;
				}
				filteredNormal += weight * face_area[currentIdx] * previous_normals[currentIdx];
			}
			if (!faceNeighbor.empty())
				filtered_normals[index] = filteredNormal.normalize();
		}

		// immediate update vertex position
		updateVertexPosition(mesh, filtered_normals, vertex_iteration_number, false);

		////
		//if (iter % 4 == 0)
		//{
		//	data_manager_->setMesh(mesh);
		//	std::stringstream ss;  ss << iter << ".obj";
		//	data_manager_->ExportMeshToFile(ss.str());
		//}
	}
}