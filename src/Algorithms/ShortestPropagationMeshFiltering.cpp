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

