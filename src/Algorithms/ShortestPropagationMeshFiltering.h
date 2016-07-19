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

	void getVertexBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getRadiusBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, double radius, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getAllFaceNeighborGMNF(TriMesh &mesh, FaceNeighborType face_neighbor_type, double radius, bool include_central_face, std::vector< std::vector<TriMesh::FaceHandle> > &all_face_neighbor);

	double GaussianWeight(double distance, double sigma);
	double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);


	double getRadius(double multiple, TriMesh &mesh);
	//得到SigmaS 和SigmaR(//SigmaR指的是法向差的范围，这个没有给出某种统计值，指导的这篇默认用的是0.35
	//)
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

	void BuildGraph3(TriMesh &mesh, FaceNeighborType face_neighbor_type, std::vector<TriMesh::Normal> &normals, 
		std::vector<std::vector<graphPair> > &fromMesh);
	void BuildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, std::vector<std::vector<graphPair> > &fromMesh, std::vector<int> &swapnumber);
	//void dijkstraPath(const std::vector< std::vector<dualEdge> > &graph,
		//TriMesh::FaceHandle source, TriMesh::FaceHandle target, std::vector<TriMesh::FaceHandle> &face_path);
	void dijkstraPath3(std::vector<std::vector<graphPair> > &fromMesh, Dijkstra mesh_graph,
		int source, int target, std::vector<int> &face_path);


	void choosePath(TriMesh &mesh, FaceNeighborType face_neighbor_type,
		TriMesh::FaceHandle fh_start, TriMesh::FaceHandle fh_end, TriMesh::Normal direct_center,
		std::vector<TriMesh::FaceHandle> &face_path);

	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);

};



#endif