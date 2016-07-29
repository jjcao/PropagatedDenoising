#ifndef MESHDENOISINGBASE_H
#define MESHDENOISINGBASE_H

#include "../datamanager.h"
#include "../parameterset.h"

class MeshDenoisingBase
{
public:
	enum FaceNeighborType { kVertexBased, kEdgeBased, kRadiusBased, kFaceRingBased };
	enum DenoiseType { kLocal, kGlobal };
public:
    MeshDenoisingBase(DataManager *_data_manager, ParameterSet *_parameter_set);
    ~MeshDenoisingBase() {}

public:
    virtual void denoise() = 0;
    virtual void initParameters() = 0;
	
	/////////////////////////////////////////////////////////////////
	// error metric
	/////////////////////////////////////////////////////////////////
	double getMeanSquareAngleError(TriMesh &DenoisedMesh, TriMesh &OriginalMesh);

protected:
	virtual void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals) {}
	void updateVertexPosition(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary = false);


	/////////////////////////////////////////////////////////////////
    // some basic function
	/////////////////////////////////////////////////////////////////
    double getAverageEdgeLength(TriMesh &mesh);
	//以SigmaS(面心距离平均值)为基准来决定取多大的半径(也即是能够选多少的邻居)
	//默认取得是2倍的面心平均距离
	double getAveragefaceCenterDistance(TriMesh &mesh);
    void getFaceArea(TriMesh &mesh, std::vector<double> &area);
    void getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid);
	void getFaceNormal(TriMesh &mesh, std::vector<TriMesh::Normal> &normals);
	void checkBadFace(TriMesh &mesh);

	/////////////////////////////////////////////////////////////////
    // two stage method, first update normals, second update vertices
	/////////////////////////////////////////////////////////////////    
    void getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, FaceNeighborType face_neighbor_type, 
		std::vector<TriMesh::FaceHandle> &face_neighbor, double radius = 0.0);
	void getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, FaceNeighborType face_neighbor_type,
		std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, double radius = 0.0);

    void getAllFaceNeighbor(TriMesh &mesh, std::vector< std::vector<TriMesh::FaceHandle> > &all_face_neighbor, 
		FaceNeighborType face_neighbor_type = kVertexBased, bool include_central_face = false, double radius=0.0);
	void getAllFaceNeighbor(TriMesh &mesh, std::vector< std::vector<std::vector<TriMesh::FaceHandle> > > &all_face_neighbor,
		FaceNeighborType face_neighbor_type = kFaceRingBased, bool include_central_face = false, double radius = 0.0);

private:
	void getVertexBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getRadiusBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, double radius, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getFaceRingBasedFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, double radius, std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor);

public:
    DataManager *data_manager_;
    ParameterSet *parameter_set_;
};

#endif // MESHDENOISINGBASE_H
