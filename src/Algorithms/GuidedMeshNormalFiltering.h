#ifndef GUIDEDMESHNORMALFILTERING_H
#define GUIDEDMESHNORMALFILTERING_H

/*
 @brief: "Guided Mesh Normal Filtering" class
 @reference: Guided Mesh Normal Filtering, PG2015
*/

#include "MeshDenoisingBase.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <queue>
#include <utility>

class GuidedMeshNormalFiltering : public MeshDenoisingBase
{
public:
    GuidedMeshNormalFiltering(DataManager *_data_manager, ParameterSet *_parameter_set);
    ~GuidedMeshNormalFiltering() {}

private:
    void denoise();
    void initParameters();


    double GaussianWeight(double distance, double sigma);
    double NormalDistance(const TriMesh::Normal &n1, const TriMesh::Normal &n2);

    void getFaceNeighborInnerEdge(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &face_neighbor, std::vector<TriMesh::EdgeHandle> &inner_edge);
    void getRangeAndMeanNormal(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                                std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
                                std::vector<std::pair<double, TriMesh::Normal> > &range_and_mean_normal);
    void getGuidedNormals(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_guided_neighbor,
                          std::vector<double> &face_area, std::vector<TriMesh::Normal> &normals,
                          std::vector<std::pair<double, TriMesh::Normal> > range_and_mean_normal,
                          std::vector<TriMesh::Normal> &guided_normals);

    double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

    void updateFilteredNormals(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
    void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
    void updateFilteredNormalsGlobalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);
};

#endif // GUIDEDMESHNORMALFILTERING_H
