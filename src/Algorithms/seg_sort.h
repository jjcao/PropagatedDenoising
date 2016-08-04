#ifndef SEG_SORT_H
#define SEG_SORT_H

#include "MeshDenoisingBase.h"

struct Projection
{
	TriMesh::Point p_point;     //projective point
	int face_index;             //original point index;
	double dist;                //geometry distance to center point
	int ring;                   //ring(topology) distance to center point
};

struct Pathmark
{
	int target;// initialized -1;
	std::vector<int> path;
};


class Segmen
{
public:
	//---------------------------FaceRing interface----------------------//
	/*the first two parameters define the original points, the third parameter define the projectoin normal, the last is used to calculating the dist*/
	void localProjection(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, const std::vector<TriMesh::Point> &face_centroid, 
		const std::vector<TriMesh::Normal> &face_normals, int centerf_index);

	//---------------------------Radius interface----------------------//
	void localProjection(const std::vector<TriMesh::FaceHandle> &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
		const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid);

	/*three points segment the place to seven regions*/
	void segmentation(const std::vector<TriMesh::Point> &fpoint);

	/*sort*/
	void calculateSort();

	/* calculate all paths*/
	void calculateAllPath(int centerf_index, std::vector<Pathmark> &pathms);

	// calculate adaptive sigma
	double calculateAdaptiveSigma(const std::vector<TriMesh::Normal> &face_normals, double smoothness);
	void clearSegmen()
	{
		projection_sets.clear();
		regions.clear();
	}

private:
	std::vector<Projection> projection_sets;
	std::vector< std::vector<Projection> > regions;

	///* ascend */
	//bool Comp1(const Projection &proc1, const Projection &proc2)
	//{
	//	if (proc1.ring == proc2.ring)
	//	{
	//		return proc1.dist < proc2.dist;
	//	}
	//	else
	//		return proc1.ring < proc2.ring;
	//}
	//bool Comp2(const Projection &proc1, const Projection &proc2)
	//{
	//	return proc1.ring < proc2.ring;
	//}
	//bool Comp3(const Projection &proc1, const Projection &proc2)
	//{
	//	return proc1.dist < proc2.dist;
	//}

	void calculatePath(std::vector<Projection> &region, int tail, int centerf_index, std::vector<Pathmark> &paths, Pathmark &pathm);
	void midPath(int mid, int centerf_index, std::vector<Projection> &region, std::vector<int> &path);


};


#endif