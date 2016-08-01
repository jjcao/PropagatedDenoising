#ifndef PROJECTIVEPROPAGATIONMESHFILTERING_H
#define PROJECTIVEPROPAGATIONMESHFILTERING_H

#include "PropagationMeshFiltering.h"
#include "seg_sort.h"

class ProjectivePropagationMeshFiltering : public PropagationMeshFiltering
{
public:
	ProjectivePropagationMeshFiltering(DataManager *_data_manager, ParameterSet *_parameter_set) : PropagationMeshFiltering(_data_manager, _parameter_set){}
	~ProjectivePropagationMeshFiltering(){}

private:
	/*the first two parameters define the original points, the third parameter define the projectoin normal, the last is used to calculating the dist*/
	void calculateGlobalPath(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, 
		const std::vector<TriMesh::Point> &face_centroid, const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid, int centerf_index, 
		const std::vector<TriMesh::Point> &fpoint, std::vector<Pathmark> &pathms);

	virtual double calculateSigma(const std::vector<TriMesh::Normal> &face_normals, TriMesh::FaceIter sourceFaceIter, double smoothness);
	
	virtual void setAllFaceNeighbor(TriMesh &mesh, FaceNeighborType face_neighbor_type, bool include_central_face, double radius);
	// facePaths[i] is the ith path from some target to source
	virtual void computeGlobalPath(TriMesh &mesh, TriMesh::FaceIter sourceFaceIter,
		const std::vector<TriMesh::Point>& face_centroid, const std::vector<TriMesh::Normal>& face_normals, std::vector<std::vector<int> > &facePaths);

	Segmen _seg;
	std::vector< std::vector<std::vector<TriMesh::FaceHandle> > > _allFaceNeighbor;
};

#endif
