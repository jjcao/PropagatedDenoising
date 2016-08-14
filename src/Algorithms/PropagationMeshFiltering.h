#pragma once
#include "MeshDenoisingBase.h"
class PropagationMeshFiltering :
	public MeshDenoisingBase
{
public:
	PropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set);
	~PropagationMeshFiltering(){}
	void denoise();
	virtual void initParameters();

protected:
	//自适应局部区域的法线情况，计算法线的方差，用方差来作为高斯的方差，（注高斯函数中的sigma是标准差）
	virtual double calculateSigma(const std::vector<TriMesh::Normal> &face_normals, TriMesh::FaceIter sourceFaceIter, int iter, double smoothness) = 0;
	virtual void setAllFaceNeighbor(TriMesh &mesh, FaceNeighborType face_neighbor_type, bool include_central_face, double radius) = 0;	
	// facePaths[i] is the ith path from some target to source


	void getFaceNeighborInnerEdge(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &face_neighbor, std::vector<TriMesh::EdgeHandle> &inner_edge);
	void getRangeAndMeanNormal(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_guided_neighbor,
		std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
		std::vector<std::pair<double, TriMesh::Normal> > &range_and_mean_normal);
	virtual void getGuidedNormals(TriMesh &mesh,
		std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
		std::vector<std::pair<double, TriMesh::Normal> > range_and_mean_normal, std::vector<TriMesh::Normal> &guided_normals) = 0;


	virtual void computeGlobalPath(TriMesh &mesh, TriMesh::FaceIter sourceFaceIter,
		const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals, std::vector<std::vector<int> > &facePaths) = 0;

	virtual void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);

	std::vector<int> _local2GlobalIdx;
	std::map<int, int> _global2localIdx;
	double GaussianWeight(double distance, double sigma);
	double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);

private:
	void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
	//void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);

	//double GaussianWeight(double distance, double sigma);
	//double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);
	
	//网格上所有面之间的面心距离的平均值的倍数，作为衡量基准，	//一般取一倍，所以值就是面心距离的平均值
	double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);
	double calculateSigmaSandR(std::vector<double> &datas, double smoothness);
	//double calculateSigma(std::vector<TriMesh::Normal> &face_normals, std::vector<TriMesh::FaceHandle> &faceNeighbor, double smoothness);
};

