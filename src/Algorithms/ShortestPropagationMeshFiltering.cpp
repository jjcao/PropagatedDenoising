#include"ShortestPropagationMeshFiltering.h"
#include <numeric>

void ShortestPropagationMeshFiltering::initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, 
	int sourceIdxGlobal, const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals)
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
		TriMesh::Point c1 = face_centroid[iter->idx()];// mesh.calc_face_centroid(*iter);
		TriMesh::Normal n1 = face_normals[iter->idx()];
		std::vector<TriMesh::FaceHandle> locFaceNeighbor;
		// todo edgeBased v.s. vertexBased
		getFaceNeighbor(mesh, *iter, kEdgeBased, locFaceNeighbor);

		for (std::vector<TriMesh::FaceHandle>::iterator it = locFaceNeighbor.begin(); it != locFaceNeighbor.end(); ++it)
		{
			TriMesh::Point c2 = face_centroid[it->idx()];
			TriMesh::Normal n2 = face_normals[it->idx()];
			// todo how about use normal distance as edge weights.
			//_localGraph[i].push_back(std::make_pair(_global2localIdx[it->idx()], std::max((n1 - n2).length(), std::numeric_limits<double>::epsilon()) ) );
			_localGraph[i].push_back(std::make_pair(_global2localIdx[it->idx()], (c1 - c2).length()));
		}
	}

	std::vector<double> distance;
	_dij.computeDistances(_localGraph, _global2localIdx[sourceIdxGlobal], distance);
}

double ShortestPropagationMeshFiltering::calculateSigma(const std::vector<TriMesh::Normal> &face_normals, 
	TriMesh::FaceIter sourceFaceIter, int iter, double smoothness)
{
	double sigma(0.25*1.4142*smoothness); // make Gaussian for angle > pi/2 == 0.0;
	sigma = std::max(0.1, sigma - 0.015*iter*2);
	return sigma;


	//////////////////////////////
	int sourceIdxGlobal = sourceFaceIter->idx();
	std::vector<TriMesh::FaceHandle> &faceNeighbor = _allFaceNeighbor[sourceIdxGlobal];

	TriMesh::Normal aver_local_normal(0.0, 0.0, 0.0);
	int len = faceNeighbor.size();
	for (int st = 0; st < len; st++)
	{
		aver_local_normal += face_normals[faceNeighbor[st].idx()];
	}
	aver_local_normal.normalize();

	double stdard = 0.0;
	for (int st = 0; st < len; st++)
	{
		double dtemp = (aver_local_normal - face_normals[faceNeighbor[st].idx()]).length();
		stdard += dtemp * dtemp;
		//stdard += NormalDistance(aver_local_normal, previous_normals[projections[i][st].face_index]);
	}

	sigma = sqrt(stdard / len)+0.01;//sigma_s替代， 作为光滑的参数; + smoothness

	return sigma;
}
void ShortestPropagationMeshFiltering::setAllFaceNeighbor(TriMesh &mesh, FaceNeighborType face_neighbor_type, bool include_central_face, double radius)
{
	getAllFaceNeighbor(mesh, _allFaceNeighbor, face_neighbor_type, include_central_face, radius);
}
void ShortestPropagationMeshFiltering::computeGlobalPath(TriMesh &mesh, TriMesh::FaceIter sourceFaceIter,
	const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals, std::vector<std::vector<int> > &facePaths)
{
	int sourceIdxGlobal = sourceFaceIter->idx();
	std::vector<TriMesh::FaceHandle> &faceNeighbor = _allFaceNeighbor[sourceIdxGlobal];
	initLocalGraph(mesh, faceNeighbor, sourceIdxGlobal, face_centroid, face_normals);

	facePaths.resize(faceNeighbor.size());
	for (int j = 0; j < (int)faceNeighbor.size(); ++j)
	{
		int targetIdxGlobal = faceNeighbor[j].idx();
		std::vector<int>& facePath = facePaths[j];
		_dij.computePath(_localGraph, _global2localIdx[targetIdxGlobal], facePath);
		for (int k = 0; k < facePath.size(); ++k)
			facePath[k] = _local2GlobalIdx[facePath[k]];
	}
}
