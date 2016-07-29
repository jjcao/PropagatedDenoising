#ifndef SHORTESTPROPAGATIONMESHFILTERING_H
#define SHORTESTPROPAGATIONMESHFILTERING_H

#include "PropagationMeshFiltering.h"
#include <map>
#include "dijkstra.h"

class ShortestPropagationMeshFiltering : public PropagationMeshFiltering
{
public:
	ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set) : PropagationMeshFiltering(_data_manager, _parameter_set){}
	~ShortestPropagationMeshFiltering(){}	

private:
	virtual void initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, int sourceIdxGlobal);
	virtual void computeGlobalPath(int targetIdxGlobal, std::vector<int>& facePath);

	virtual void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);

	Dijkstra _dij;
	AdjacencyList _localGraph;	

};



#endif