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

	//�õ�SigmaS ��SigmaR(//SigmaRָ���Ƿ����ķ�Χ�����û�и���ĳ��ͳ��ֵ��ָ������ƪĬ���õ���0.35
	//)
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

	void buildLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> face_neighbor, AdjacencyList &localGraph, 
		std::vector<int>& local2GlobalIdx, std::map<int, int>& global2localIdx);

	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
};



#endif