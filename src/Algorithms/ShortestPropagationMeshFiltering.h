#ifndef SHORTESTPROPAGATIONMESHFILTERING_H
#define SHORTESTPROPAGATIONMESHFILTERING_H

#include "MeshDenoisingBase.h"
#include "dijkstra.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <queue>

#include <limits>
#include <algorithm>
#include <iostream>


class ShortestPropagationMeshFiltering : public MeshDenoisingBase
{
public:
	ShortestPropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set);
	~ShortestPropagationMeshFiltering(){}

private:
	void denoise();
	void initParameters();

	double GaussianWeight(double distance, double sigma);
	double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);

	//得到SigmaS 和SigmaR(//SigmaR指的是法向差的范围，这个没有给出某种统计值，指导的这篇默认用的是0.35
	//)
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

	void buildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, AdjacencyList &localGraph, std::vector<int> &localFaceIdx);
	//void dijkstraPath(const std::vector< std::vector<dualEdge> > &graph,
		//TriMesh::FaceHandle source, TriMesh::FaceHandle target, std::vector<TriMesh::FaceHandle> &face_path);
	void dijkstraPath3(std::vector<std::vector<GraphPair> > &fromMesh, Dijkstra mesh_graph,
		int source, int target, std::vector<int> &face_path);


	void choosePath(TriMesh &mesh, FaceNeighborType face_neighbor_type,
		TriMesh::FaceHandle fh_start, TriMesh::FaceHandle fh_end, TriMesh::Normal direct_center,
		std::vector<TriMesh::FaceHandle> &face_path);

	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
};



#endif