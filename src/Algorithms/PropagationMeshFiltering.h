#pragma once
#include "MeshDenoisingBase.h"
class PropagationMeshFiltering :
	public MeshDenoisingBase
{
public:
	PropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set);
	~PropagationMeshFiltering(){}
	void denoise();

protected:
	virtual void initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, int sourceIdxGlobal) = 0;
	virtual void computeGlobalPath(int targetIdxGlobal, std::vector<int>& facePath) = 0;

	virtual void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals) = 0;

	std::vector<int> _local2GlobalIdx;
	std::map<int, int> _global2localIdx;
	double GaussianWeight(double distance, double sigma);

private:
	void initParameters();

	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	//void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);

	//double GaussianWeight(double distance, double sigma);
	double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);
	
	//网格上所有面之间的面心距离的平均值的倍数，作为衡量基准，	//一般取一倍，所以值就是面心距离的平均值
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);	
	//double calculateSigma(std::vector<TriMesh::Normal> &face_normals, std::vector<TriMesh::FaceHandle> &faceNeighbor, double smoothness);
};

