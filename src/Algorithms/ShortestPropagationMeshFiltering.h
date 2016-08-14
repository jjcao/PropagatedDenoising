#ifndef SHORTESTPROPAGATIONMESHFILTERING_H
#define SHORTESTPROPAGATIONMESHFILTERING_H

#include "PropagationMeshFiltering.h"
#include "dijkstra.h"

class ShortestPropagationMeshFiltering : public PropagationMeshFiltering
{
public:
	ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set) : PropagationMeshFiltering(_data_manager, _parameter_set){}
	~ShortestPropagationMeshFiltering(){}	

private:
	virtual void setAllFaceNeighbor(TriMesh &mesh, FaceNeighborType face_neighbor_type, bool include_central_face, double radius);

	virtual void getGuidedNormals(TriMesh &mesh,
		std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
		std::vector<std::pair<double, TriMesh::Normal> > range_and_mean_normal, std::vector<TriMesh::Normal> &guided_normals);

	virtual void computeGlobalPath(TriMesh &mesh, TriMesh::FaceIter sourceFaceIter,
		const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals, std::vector<std::vector<int> > &facePaths);


	void initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, int sourceIdxGlobal, const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals);
	virtual void checkBadFace(TriMesh &mesh) {};
	virtual double calculateSigma(const std::vector<TriMesh::Normal> &face_normals, TriMesh::FaceIter sourceFaceIter, int iter, double smoothness);

	Dijkstra _dij;
	AdjacencyList _localGraph;	
	std::vector<std::vector<TriMesh::FaceHandle> > _allFaceNeighbor;
	std::vector<std::vector<TriMesh::FaceHandle> > _allGuidedNeighbor;


};



#endif