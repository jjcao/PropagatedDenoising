#include "ProjectivePropagationMeshFiltering.h"
void ProjectivePropagationMeshFiltering::initParameters()
{
	PropagationMeshFiltering::initParameters();
	parameter_set_->setStringListIndex("Face Neighbor", 3);
}

void ProjectivePropagationMeshFiltering::setAllFaceNeighbor(TriMesh &mesh, FaceNeighborType face_neighbor_type, bool include_central_face, double radius)
{
	getAllFaceNeighbor(mesh, _allFaceNeighbor, face_neighbor_type, include_central_face, radius);	
}

void ProjectivePropagationMeshFiltering::computeGlobalPath(TriMesh &mesh, TriMesh::FaceIter sourceFaceIter, 
	const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal> &face_normals, std::vector<std::vector<int> > &facePaths)
{
	//存储当前面片的三个顶点坐标，为分区域做好准备
	std::vector<TriMesh::Point> fpoint;
	fpoint.resize(3); int pindex = 0;
	for (TriMesh::FaceVertexIter fv_it = mesh.fv_iter(*sourceFaceIter); fv_it.is_valid(); fv_it++)
	{
		fpoint[pindex] = mesh.point(*fv_it);
		pindex++;
	}

	int centerf_index = sourceFaceIter->idx();

	std::vector<Pathmark> pathms;
	//--------FaceRingBased------//
	std::vector<std::vector<TriMesh::FaceHandle> > &face_neighbor = _allFaceNeighbor[centerf_index];
	calculateGlobalPath(face_neighbor, face_centroid, face_normals, centerf_index, fpoint, pathms);

	facePaths.reserve(pathms.size());
	for (int i = 0; i < pathms.size(); ++i)
		//facePaths[i] = pathms[i].path;
		facePaths.push_back(pathms[i].path);
}

void ProjectivePropagationMeshFiltering::calculateGlobalPath(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, 
	const std::vector<TriMesh::Point> &face_centroid, const std::vector<TriMesh::Normal> &face_normals, int centerf_index,
	const std::vector<TriMesh::Point> &fpoint, std::vector<Pathmark> &pathms)
{
	_seg.localProjection(face_neighbor, face_centroid, face_normals, centerf_index);
	_seg.segmentation(fpoint);
	_seg.calculateSort();
	_seg.calculateAllPath(centerf_index, pathms);
}

double ProjectivePropagationMeshFiltering::calculateSigma(const std::vector<TriMesh::Normal> &face_normals, 
	TriMesh::FaceIter sourceFaceIter, int iter, double smoothness)
{
	double sigma = _seg.calculateAdaptiveSigma(face_normals, smoothness);
	_seg.clearSegmen();

	return sigma;
}
