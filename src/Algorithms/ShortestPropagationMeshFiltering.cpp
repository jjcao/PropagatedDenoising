#include"ShortestPropagationMeshFiltering.h"
#include<iostream>


ShortestPropagationMeshFiltering::ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set)
	: MeshDenoisingBase(_data_manager, _parameter_set)
{
	initParameters();
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

//  如果使用则使用face_neighbor 来选择使用哪一种邻居
void ShortestPropagationMeshFiltering::getVertexBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, std::vector<TriMesh::FaceHandle> &face_neighbor)
{
	//getFaceNeighbor(mesh, fh, kVertexBased, face_neighbor);
	getFaceNeighbor(mesh, fh, kEdgeBased, face_neighbor);
}

//得到面中心的所有小于给定半径值的面片
void ShortestPropagationMeshFiltering::getRadiusBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, double radius, std::vector<TriMesh::FaceHandle> &face_neighbor)
{
	TriMesh::Point ci = mesh.calc_face_centroid(fh);
	std::vector<bool> flag((int)mesh.n_faces(), false);

	face_neighbor.clear();
	flag[fh.idx()] = true;
	std::queue<TriMesh::FaceHandle> queue_face_handle;
	queue_face_handle.push(fh);

	//这种迭代策略一定要熟悉
	std::vector<TriMesh::FaceHandle> temp_face_neighbor;
	while (!queue_face_handle.empty())
	{
		TriMesh::FaceHandle temp_face_handle_queue = queue_face_handle.front();
		if (temp_face_handle_queue != fh)
			face_neighbor.push_back(temp_face_handle_queue);
		queue_face_handle.pop();
		getVertexBasedFaceNeighbor(mesh, temp_face_handle_queue, temp_face_neighbor);
		for (int i = 0; i < (int)temp_face_neighbor.size(); i++)
		{
			TriMesh::FaceHandle temp_face_handle = temp_face_neighbor[i];
			if (!flag[temp_face_handle.idx()])
			{
				TriMesh::Point cj = mesh.calc_face_centroid(temp_face_handle);
				double distance = (ci - cj).length();
				if (distance <= radius)
					queue_face_handle.push(temp_face_handle);
				flag[temp_face_handle.idx()] = true;
			}
		}
	}
}

void ShortestPropagationMeshFiltering::getAllFaceNeighborGMNF(TriMesh &mesh, MeshDenoisingBase::FaceNeighborType face_neighbor_type, double radius, bool include_central_face,
	std::vector<std::vector<TriMesh::FaceHandle> > &all_face_neighbor)
{
	std::vector<TriMesh::FaceHandle> face_neighbor;
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		if (face_neighbor_type == kVertexBased)//--   1-ring
			getVertexBasedFaceNeighbor(mesh, *f_it, face_neighbor);
		else if (face_neighbor_type == kRadiusBased)//--  r-ring
			getRadiusBasedFaceNeighbor(mesh, *f_it, radius, face_neighbor);

		if (include_central_face)
			face_neighbor.push_back(*f_it);
		all_face_neighbor[f_it->idx()] = face_neighbor;
	}
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

//以SigmaS(面心距离平均值)为基准来决定取多大的半径(也即是能够选多少的邻居)
//默认取得是2倍的面心平均距离
double ShortestPropagationMeshFiltering::getRadius(double multiple, TriMesh &mesh)
{
	std::vector<TriMesh::Point> centroid;
	getFaceCentroid(mesh, centroid);

	double radius = 0.0;
	double num = 0.0;
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point fi = centroid[f_it->idx()];
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
		{
			TriMesh::Point fj = centroid[ff_it->idx()];
			radius += (fj - fi).length();//    面中心距离的平均值
			num++;
		}
	}
	return radius * multiple / num;
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

void ShortestPropagationMeshFiltering::BuildGraph3(TriMesh &mesh, FaceNeighborType face_neighbor_type, std::vector<TriMesh::Normal> &normals,
	std::vector<std::vector<graphPair> > &fromMesh)
{
	fromMesh.clear();
	getFaceNormal(mesh, normals);
	fromMesh.resize(mesh.n_faces());
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		//std::vector<graphPair> G_edges;
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_iter(f_it); ff_it.is_valid(); ff_it++)
		{
			//graph[id1].push_back(std::make_pair(id2, weight));
			//fromMesh[f_it->idx()].push_back(std::vector<graphPair>::value_type(ff_it->idx(), NormalDistance(normals[(*f_it).idx()], normals[(*ff_it).idx()])));
			fromMesh[f_it->idx()].push_back(std::make_pair(ff_it->idx(), NormalDistance(normals[(*f_it).idx()], normals[(*ff_it).idx()])));
		}
	}

}

int match_number(std::vector<int> &swapnumber, int value)
{
	int i = 0;
	for (std::vector<int>::iterator iter = swapnumber.begin(); iter != swapnumber.end(); iter++)
	{
		if (*iter == value)
			return i;
		else
			i++;
	}
	return -1;
}
void match_face_path(std::vector<int> &swapnumber, std::vector<int> &face_path, std::vector<int> &match_facep)
{
	match_facep.clear();
	for (std::vector<int>::iterator iter = face_path.begin(); iter != face_path.end(); iter++)
	{
		match_facep.push_back(swapnumber[*iter]);
	}
}
void ShortestPropagationMeshFiltering::BuildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, std::vector<std::vector<graphPair> > &fromMesh, std::vector<int> &swapnumber)
{
	fromMesh.clear(); swapnumber.clear();
	fromMesh.resize(face_neighbor.size());
	//std::vector<int> swapnumber;
	for (std::vector<TriMesh::FaceHandle>::iterator iter = face_neighbor.begin(); iter != face_neighbor.end(); iter++)
	{
		swapnumber.push_back(iter->idx());
	}
	sort(swapnumber.begin(),swapnumber.end());
	auto end_unique = unique(swapnumber.begin(), swapnumber.end());
	swapnumber.erase(end_unique, swapnumber.end());
	for (std::vector<TriMesh::FaceHandle>::iterator iter = face_neighbor.begin(); iter != face_neighbor.end(); iter++)
	{
		TriMesh::Point c1 = mesh.calc_face_centroid(*iter);
		int i = match_number(swapnumber, iter->idx());

		for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(*iter); fv_it.is_valid(); fv_it++)
		{
			for (TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*fv_it); vf_it.is_valid(); vf_it++)
			{
				TriMesh::Point c2 = mesh.calc_face_centroid(*vf_it);
				int j = match_number(swapnumber, vf_it->idx());
				if (j == -1)
				{
					break;
				}
				fromMesh[i].push_back(std::make_pair(j, (c1-c2).length()));
			}
		}
	}

}
void ShortestPropagationMeshFiltering::dijkstraPath3(std::vector<std::vector<graphPair> > &fromMesh, Dijkstra mesh_graph, int source, int target, std::vector<int> &face_path)
{
	face_path.clear();
	//Dijkstra mesh_graph;
	mesh_graph.buildGraph(fromMesh);
	mesh_graph.initialization(source);
	mesh_graph.computeDistances(target);
	mesh_graph.createPath(source, target);
	//face_path = mesh_graph.path;
	for (std::vector<int>::iterator iter = mesh_graph.path.begin(); iter != mesh_graph.path.end(); iter++)
	{
		face_path.push_back(*iter);
	}
}


void ShortestPropagationMeshFiltering::choosePath(TriMesh &mesh, FaceNeighborType face_neighbor_type,
	TriMesh::FaceHandle fh_start, TriMesh::FaceHandle fh_end, TriMesh::Normal direct_center,
	std::vector<TriMesh::FaceHandle> &face_path)
{
	face_path.clear();
	face_path.push_back(fh_start);
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	getFaceCentroid(mesh, face_centroid);
	int index_end = fh_end.idx();
	bool bad_path = false;

	while (fh_start != fh_end)
	{
		std::vector<TriMesh::FaceHandle> face_neighbor_cycle;
		face_neighbor_cycle.clear();
		//getFaceNeighbor(mesh, fh_start, face_neighbor_type, face_neighbor_cycle);
		getVertexBasedFaceNeighbor(mesh, fh_start, face_neighbor_cycle);

		//find the minimum angle have some question, emerging exceed the goal
		//apply the minimum distance

		int index_start = fh_start.idx();
		double minDistance = 1e6;

		TriMesh::FaceHandle mid_fh;
		for (int i = 0; i < (int)face_neighbor_cycle.size(); i++)
		{
			int index_mid = face_neighbor_cycle[i].idx();
			double value_mid = (face_centroid[index_mid] - face_centroid[index_end]).length();
			if (value_mid < minDistance)
			{
				minDistance = value_mid;
				mid_fh = face_neighbor_cycle[i];
			}

		}
		face_path.push_back(mid_fh);
		fh_start = mid_fh;
		if (face_path.size() > 10)
			break;

	}

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
	double multiple_sigma_s;
	if (!parameter_set_->getValue(QString("Multiple(* sigma_s)"), multiple_sigma_s))
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

	//FaceNeighborType face_neighbor_type = face_neighbor_index == 0 ? kRadiusBased : kVertexBased;
	FaceNeighborType face_neighbor_type = kRadiusBased;
	double radius;
	if (face_neighbor_type == kRadiusBased)
		radius = getRadius(multiple_radius, mesh);

	std::vector<std::vector<TriMesh::FaceHandle> > all_face_neighbor((int)mesh.n_faces());
	getAllFaceNeighborGMNF(mesh, face_neighbor_type, radius, include_central_face, all_face_neighbor);

	getFaceNormal(mesh, filtered_normals);

	std::vector<double> face_area((int)mesh.n_faces());
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	std::vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());//上一次的法线

	std::ofstream fout("face_path.txt");
	for (int iter = 0; iter < normal_iteration_number; iter++)
	{
		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		getFaceNormal(mesh, previous_normals);

		//---------------------------          构建全局图          -----------------------------
		//std::vector<std::vector<graphPair> > fromGlobalMesh;
		//BuildGraph3(mesh, kVertexBased, previous_normals, fromGlobalMesh);//基于边的邻居，应该改为基于顶点的
		//Dijkstra mesh_graph;
		//----------------------------------------------------------------------------------

		//get SigamS
		//get SigmaR
		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			int index = f_it->idx();
			std::vector<TriMesh::FaceHandle> face_neighbor = all_face_neighbor[index];

			//---------------------------------------得到局部图-------------------------------------
			std::vector<std::vector<graphPair> > fromMesh;
			std::vector<int> swapnumber;
			BuildLocalGraph(mesh, face_neighbor, fromMesh, swapnumber);
			Dijkstra mesh_graph;
			//--------------------------------------------------------------------------------------

			TriMesh::Normal filtered_normal(0.0, 0.0, 0.0);
			for (int j = 0; j < (int)face_neighbor.size(); j++)
			{
				int current_face_index = face_neighbor[j].idx();
				//TriMesh::Normal direct_center = (face_centroid[index] - face_centroid[current_face_index]).normalize();
				std::vector<int> face_path;
				double weight = 0.0;
				if (face_neighbor[j] != *f_it)
				{
					//choosePath(mesh, face_neighbor_type, face_neighbor[j], *f_it, direct_center, face_path);

					//dijkstraPath3(fromMesh, mesh_graph, (*f_it).idx(), (face_neighbor[j]).idx(), face_path);

					int source = match_number(swapnumber, (*f_it).idx());
					int target = match_number(swapnumber, (face_neighbor[j]).idx());
					dijkstraPath3(fromMesh, mesh_graph, source, target, face_path);
					std::vector<int> match_facep;
					match_face_path(swapnumber, face_path, match_facep);
					face_path.swap(match_facep);
					for (std::vector<int>::iterator iterpp = face_path.begin(); iterpp != face_path.end(); iterpp++)
					{
						fout << (*iterpp) << " ";
					}
					fout << std::endl;

					//fout << face_path.size() << std::endl;

					double sumPF1 = 0.0, sumPF2 = 0.0;
					//前后两项的法线差异
					for (int k = (int)face_path.size() - 1; k > 0; k--)
					{
						int preIndex = face_path[k];
						int nexIndex = face_path[k - 1];
						double temp11 = (previous_normals[nexIndex] - previous_normals[preIndex]).length();
						sumPF1 += temp11*temp11;
					}
					//累积路径差异
					for (int k = (int)face_path.size() - 1; k >= 0; k--)
					{
						int currentIndex = face_path[k];
						double temp22 = (previous_normals[currentIndex] - previous_normals[index]).length();
						sumPF2 += temp22*temp22;
					}

					//double sigmaA = 1, sigmaR = 1;
					weight = GaussianWeight(sqrt(sumPF1), multiple_sigma_s) * GaussianWeight(sqrt(sumPF2), sigma_r);
				}
				else
				{
					weight = 1.0;
				}
				filtered_normal += weight * face_area[current_face_index] * previous_normals[current_face_index];
			}
			if (face_neighbor.size())
				filtered_normals[index] = filtered_normal.normalize();
		}
		// immediate update vertex position
		updateVertexPosition(mesh, filtered_normals, vertex_iteration_number, false);
	}
	fout.close();
}
