#ifndef SHORTESTPROPAGATIONMESHFILTERING_H
#define SHORTESTPROPAGATIONMESHFILTERING_H

#include "MeshDenoisingBase.h"
#include <map>
#include "dijkstra.h"

class ShortestPropagationMeshFiltering : public MeshDenoisingBase
{
public:
	ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set);
	~ShortestPropagationMeshFiltering(){}

	void denoise();
private:
	void initParameters();

	double GaussianWeight(double distance, double sigma);
	double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);

	//得到SigmaS 和SigmaR(//SigmaR指的是法向差的范围，这个没有给出某种统计值，指导的这篇默认用的是0.35
	//)
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

	void buildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, AdjacencyList &localGraph, 
		std::vector<int>& local2GlobalIdx, std::map<int, int>& global2localIdx);

	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
};



#endif