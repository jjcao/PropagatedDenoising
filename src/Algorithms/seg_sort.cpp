#include"seg_sort.h"
#include<algorithm>

void Segmen::localProjection(const std::vector< std::vector<TriMesh::FaceHandle> > &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
	const std::vector<TriMesh::Normal> &face_normals, int centerf_index)
{
	//这个函数的效率太低了
	projection_sets.clear();

	TriMesh::Point centerf_centroid = face_centroid[centerf_index];
	TriMesh::Normal centerf_normal = face_normals[centerf_index];
	for (int j = 1; j < (int)face_neighbor.size(); ++j)
	{
		const std::vector<TriMesh::FaceHandle> &face_neighbor_ring = face_neighbor[j];
		for (int k = 0; k < face_neighbor_ring.size(); ++k)
		{
			int neighborf_index = face_neighbor_ring[k].idx();
			TriMesh::Point neighborf_centroid = face_centroid[neighborf_index];
			TriMesh::Point neighbotf_normal = face_normals[neighborf_index];

			Projection temp;
			temp.face_index = neighborf_index;
			temp.p_point = neighborf_centroid - abs(centerf_normal | (neighborf_centroid - centerf_centroid))*centerf_normal;
			//temp.dist = (neighborf_centroid - centerf_centroid).length();
			//------ using normal distance, not work well, almost the same result
			temp.dist = (neighbotf_normal - centerf_normal).length();
			temp.ring = j;

			projection_sets.push_back(temp);
		}

	}//end projection sets
}

void Segmen::localProjection(const std::vector<TriMesh::FaceHandle> &face_neighbor, const std::vector<TriMesh::Point> &face_centroid,
	const TriMesh::Normal &centerf_normal, const TriMesh::Point &centerf_centroid)
{
	projection_sets.clear();

	for (int i = 0; i < face_neighbor.size(); ++i)
	{
		int neighborf_index = face_neighbor[i].idx();
		TriMesh::Point neighborf_centroid = face_centroid[neighborf_index];
		//TriMesh::Point neighbotf_normal = normals[neighborf_index];

		Projection temp;
		temp.face_index = neighborf_index;
		temp.p_point = neighborf_centroid - abs(centerf_normal | (neighborf_centroid - centerf_centroid))*centerf_normal;
		temp.dist = (neighborf_centroid - centerf_centroid).length();
		//------ using normal distance, not work well, almost the same result
		//temp.dist = (neighbotf_normal - centerf_normal).length();

		projection_sets.push_back(temp);
	}
}

void Segmen::segmentation(const std::vector<TriMesh::Point> &fpoint)
{
	regions.clear(); regions.resize(6);

	for (std::vector<Projection>::iterator p_iter = projection_sets.begin(); p_iter != projection_sets.end(); p_iter++)
	{
		////-----six regions
		double a1 = fpoint[1][1] * fpoint[2][2] - fpoint[1][2] * fpoint[2][1];
		double b1 = fpoint[1][0] * fpoint[2][2] - fpoint[1][2] * fpoint[2][0];
		double c1 = fpoint[1][0] * fpoint[2][1] - fpoint[1][1] * fpoint[2][0];

		double a2 = fpoint[2][1] * fpoint[0][2] - fpoint[2][2] * fpoint[0][1];
		double b2 = fpoint[2][0] * fpoint[0][2] - fpoint[2][2] * fpoint[0][0];
		double c2 = fpoint[2][0] * fpoint[0][1] - fpoint[2][1] * fpoint[0][0];

		double a3 = fpoint[0][1] * fpoint[1][2] - fpoint[0][2] * fpoint[1][1];
		double b3 = fpoint[0][0] * fpoint[1][2] - fpoint[0][2] * fpoint[1][0];
		double c3 = fpoint[0][0] * fpoint[1][1] - fpoint[0][1] * fpoint[1][0];

		double x = (p_iter->p_point)[0];
		double y = (p_iter->p_point)[1];
		double z = (p_iter->p_point)[2];

		double area_coor1 = x * a1 - y * b1 + z * c1;
		double area_coor2 = x * a2 - y * b2 + z * c2;
		double area_coor3 = x * a3 - y * b3 + z * c3;
		if (area_coor1 >= 0 && area_coor2 <= 0 && area_coor3 <= 0)
		{
			regions[0].push_back(*p_iter); continue;
		}

		if (area_coor1 >= 0 && area_coor2 >= 0 && area_coor3 <= 0)
		{
			regions[1].push_back(*p_iter); continue;
		}

		if (area_coor1 <= 0 && area_coor2 >= 0 && area_coor3 <= 0)
		{
			regions[2].push_back(*p_iter); continue;
		}

		if (area_coor1 <= 0 && area_coor2 >= 0 && area_coor3 >= 0)
		{
			regions[3].push_back(*p_iter); continue;
		}

		if (area_coor1 <= 0 && area_coor2 <= 0 && area_coor3 >= 0)
		{
			regions[4].push_back(*p_iter); continue;
		}

		if (area_coor1 >= 0 && area_coor2 <= 0 && area_coor3 >= 0)
		{
			regions[5].push_back(*p_iter); continue;
		}
		if (area_coor1 >= 0 && area_coor2 >= 0 && area_coor3 >= 0)
		{
			continue;
		}//6区结束
	}//分区结束
}

/* ascend */
bool Comp1(const Projection &proc1, const Projection &proc2)
{
	if (proc1.ring == proc2.ring)
	{
		return proc1.dist < proc2.dist;
	}
	else
		return proc1.ring < proc2.ring;
}
bool Comp2(const Projection &proc1, const Projection &proc2)
{
	return proc1.ring < proc2.ring;
}
bool Comp3(const Projection &proc1, const Projection &proc2)
{
	return proc1.dist < proc2.dist;
}


void Segmen::calculateSort()
{
	for (int i = 0; i < regions.size(); ++i)
	{
		std::sort(regions[i].begin(), regions[i].end(), Comp1);
	}
}

double Segmen::calculateAdaptiveSigma(const std::vector<TriMesh::Normal> &face_normals, double smoothness)
{
	TriMesh::Normal aver_local_normal(0.0, 0.0, 0.0);
	double minsigma = 10;
	
	for (int num = 0; num < regions.size(); ++num)
	{
		int len = regions[num].size();
		if (len > 1)
		{
			for (int st = 0; st < len; ++st)
			{
				aver_local_normal += face_normals[regions[num][st].face_index];
			}
			aver_local_normal.normalize();

			double stdard = 0.0;
			for (int st = 0; st < len; st++)
			{
				double dtemp = (aver_local_normal - face_normals[regions[num][st].face_index]).length();
				stdard += dtemp * dtemp;
			}
			stdard = sqrt(stdard / len) + smoothness;//作为光滑的参数

			if (minsigma > stdard)
				minsigma = stdard;
		}
	}
	return minsigma;
}

void Segmen::midPath(int mid, int centerf_index, std::vector<Projection> &region, std::vector<int> &path)
{
	path.clear();

	for (int i = mid; i >= 0; i--)
	{
		path.push_back(region[i * i].face_index);
	}
	path.push_back(centerf_index);
}

bool operator==(const Pathmark &a1, const Pathmark &a2)
{
	return a1.target == a2.target;
}

void Segmen::calculatePath(std::vector<Projection> &region, int tail, int centerf_index, std::vector<Pathmark> &paths, Pathmark &pathm)
{
	std::vector<int> pathtemp; 
	pathtemp.clear();

	Pathmark  findtemp;
	//findtemp.target == region[tail].face_index;

	int tail_t = tail;
	
	int k = 0;
	while (tail >= k * k)
	{
		k++;
	}

	while (1)
	{
		if (tail == (k - 1) * (k - 1))
		{
			
			findtemp.target = region[tail].face_index;
			std::vector<Pathmark>::iterator iter = find(paths.begin(), paths.end(), findtemp);
			if (iter != paths.end())
			{
				pathtemp = iter->path;
				break;
			}
			else
			{
				midPath(k - 1, centerf_index, region, pathtemp);
				break;
			}
		}
		else
		{
			//pathtemp.push_back(region[tail].face_index);

			if (tail - 1 == (k - 1) * (k - 1))
			{				
				findtemp.target = region[tail - 1].face_index;
				std::vector<Pathmark>::iterator iter = find(paths.begin(), paths.end(), findtemp);
				if (iter != paths.end())
				{
					pathtemp = iter->path;
					pathtemp.insert(pathtemp.begin(), region[tail].face_index);
					break;
				}
				else
				{
					midPath(k - 1, centerf_index, region, pathtemp);
					break;
				}
			}
			else
			{
				int _tail = tail - 2;

				findtemp.target = region[_tail].face_index;
				std::vector<Pathmark>::iterator iter = find(paths.begin(), paths.end(), findtemp);
				if (iter != paths.end())
				{
					pathtemp = iter->path;
					pathtemp.insert(pathtemp.begin(), region[tail].face_index);
					break;
				}

			}
		}
	}

	pathm.path = pathtemp;
	pathm.target = region[tail_t].face_index;

}


void Segmen::calculateAllPath(int centerf_index, std::vector<Pathmark> &paths)
{ 
	paths.clear();
	int len = regions[0].size() + regions[1].size() + regions[2].size() + regions[3].size() + regions[4].size() + regions[5].size();
	paths.resize(len);
	for (int i = 0; i < paths.size(); ++i)
	{
		paths[i].target = -1;
		paths[i].path.clear();
	}

	int number = 0;
	for (int i = 0; i < regions.size(); ++i)
	{

		for (int k = 0; k < regions[i].size(); ++k)
		{
			calculatePath(regions[i], k, centerf_index, paths, paths[number]);			
			++number;
		}
	}
	//projection_sets.clear();
	//regions.clear();
}

