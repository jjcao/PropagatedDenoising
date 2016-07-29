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
	void calculateGlobalPath(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, const std::vector<TriMesh::Point> &face_centroid, 
		const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid, int centerf_index, const std::vector<TriMesh::Point> &fpoint,
		std::vector<Pathmark> &pathms);

	/*the first two parameters define the original points, the third parameter define the projectoin normal, the last is used to calculating the dist*/
	void calculateGlobalPath(const std::vector<TriMesh::FaceHandle> &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
		const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid, int centerf_index, const std::vector<TriMesh::Point> &fpoint,
		std::vector<Pathmark> &pathms);

	virtual void updateFilteredNormalsLocalScheme(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals);


	virtual void initLocalGraph(TriMesh &mesh, std::vector<TriMesh::FaceHandle> &faceNeighbor, int sourceIdxGlobal)
	{
		return;
	}
	virtual void computeGlobalPath(int targetIdxGlobal, std::vector<int>& facePath)
	{
		return;
	}

	Segmen seg;
};

#endif
