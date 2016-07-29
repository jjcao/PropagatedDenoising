#include "ProjectivePropagationMeshFiltering.h"

void ProjectivePropagationMeshFiltering::updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals)
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
	double radius = multiple_radius;
	if (face_neighbor_type == kRadiusBased)
		radius = getAveragefaceCenterDistance(mesh) * multiple_radius;

	//--------RadiusBased-------//
	//////////////////////////////
	//////////////////////////////
	//std::vector<std::vector<TriMesh::FaceHandle> > all_face_neighbor((int)mesh.n_faces());
	//getAllFaceNeighbor(mesh, all_face_neighbor, face_neighbor_type, include_central_face, radius);

	//--------FaceRingBased------//
	///////////////////////////////
	///////////////////////////////
	std::vector< std::vector<std::vector<TriMesh::FaceHandle> > > all_face_neighbor((int)mesh.n_faces());
	getAllFaceNeighbor(mesh, all_face_neighbor, face_neighbor_type, include_central_face, radius);



	getFaceNormal(mesh, filtered_normals); //----------- be initialized?---------------
	std::vector<double> face_area((int)mesh.n_faces());
	std::vector<TriMesh::Point> face_centroid((int)mesh.n_faces());
	std::vector<TriMesh::Normal> previous_normals((int)mesh.n_faces());//上一次的法线

	//////---- record the std,
	////std::ofstream fout("graph.txt");

	for (int iter = 0; iter < normal_iteration_number; iter++)
	{
		std::cout << "iteration: " << iter << std::endl;

		getFaceCentroid(mesh, face_centroid);
		getFaceArea(mesh, face_area);
		//getFacePerimeter(mesh, face_perimeter);
		getFaceNormal(mesh, previous_normals);//------- 这里的法线是重新求出来的，不是更新后的法线？----------

		for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			TriMesh::Normal filtered_normal(0.0, 0.0, 0.0);//滤波结果初始化 

			//存储当前面片的三个顶点坐标，为分区域做好准备
			std::vector<TriMesh::Point> fpoint;
			fpoint.resize(3); int pindex = 0;
			for (TriMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				fpoint[pindex] = mesh.point(*fv_it);
				pindex++;
			}

			int centerf_index = f_it->idx();
			TriMesh::Point centerf_centroid = face_centroid[centerf_index];
			TriMesh::Normal centerf_normal = previous_normals[centerf_index];

			std::vector<Pathmark> pathms;
			//--------RadiusBased-------//
			//////////////////////////////
			//////////////////////////////
			//std::vector<TriMesh::FaceHandle> face_neighbor = all_face_neighbor[centerf_index];


			//--------FaceRingBased------//
			///////////////////////////////
			///////////////////////////////
			std::vector<std::vector<TriMesh::FaceHandle> > face_neighbor = all_face_neighbor[centerf_index];



			calculateGlobalPath(face_neighbor, face_centroid, centerf_normal, centerf_centroid, centerf_index, fpoint, pathms);
			for (int i = 0; i < pathms.size(); ++i)
			{
				std::vector<int> onepath = pathms[i].path;
				double weight = 0.0;
				double sumPF1 = 0.0, sumPF2 = 0.0;

				//前后两项差异
				for (int k = (int)onepath.size() - 1; k > 0; k--)
				{
					int preIndex = onepath[k];
					int nexIndex = onepath[k - 1];
					double temp11 = (previous_normals[nexIndex] - previous_normals[preIndex]).length();
					sumPF1 += temp11*temp11;
				}
				//累积路径差异
				int centerindex = onepath[onepath.size() - 1];
				for (int k = (int)onepath.size() - 1; k > 0; k--)//不含他自己
				{
					int tailIndex = onepath[k - 1];
					double temp22 = (previous_normals[tailIndex] - previous_normals[centerindex]).length();
					sumPF2 += temp22*temp22;
				}
				/*    change tomorrow   */

				weight = GaussianWeight(sqrt(sumPF1), sigma_s) * GaussianWeight(sqrt(sumPF2), sigma_r);
				filtered_normal += weight * face_area[onepath[0]] * previous_normals[onepath[0]];
				onepath.clear();
			}

			filtered_normals[centerf_index] = filtered_normal.normalize();

		}//end all face filtering

		// immediate update vertex position
		updateVertexPosition(mesh, filtered_normals, vertex_iteration_number, false);

		checkBadFace(mesh);//对坏的面进行微小的扰动

		////
		if (iter % 4 == 0)
		{
			data_manager_->setMesh(mesh);
			std::stringstream ss;  ss << iter << ".obj";
			data_manager_->ExportMeshToFile(ss.str());
		}
	}

	//fout.close();
}

void ProjectivePropagationMeshFiltering::calculateGlobalPath(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
	const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid, int centerf_index, const std::vector<TriMesh::Point> &fpoint,
	std::vector<Pathmark> &pathms)
{
	seg.localProjection(face_neighbor, face_centroid, centerf_normal, centerf_centroid);
	seg.segmentation(fpoint);
	seg.calculateSort();
	seg.calculateAllPath(centerf_index, pathms);
}

void ProjectivePropagationMeshFiltering::calculateGlobalPath(const std::vector<TriMesh::FaceHandle> &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
	const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid, int centerf_index, const std::vector<TriMesh::Point> &fpoint,
	std::vector<Pathmark> &pathms)
{
	seg.localProjection(face_neighbor, face_centroid, centerf_normal, centerf_centroid);
	seg.segmentation(fpoint);
	seg.calculateSort();
	seg.calculateAllPath(centerf_index, pathms);
}