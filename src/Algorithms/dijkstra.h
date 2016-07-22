#ifndef DIJKSTRA_H
#define DIJKSTRA_H

/*
	compute shortest distance and path of a graph represented by adjacency list via Dijkstra algorithm

	todo: 改成1次计算所有距离的，

	jjcao @ 2016
*/
#include <vector>

typedef std::pair<int, double> GraphPair; 
typedef std::vector<std::vector<GraphPair> > AdjacencyList;// GraphPair means <index of the node, weight between the node and ... >

class Dijkstra
{
public:
	// call it first before call computePath()
	void computeDistances(const AdjacencyList & graph, const int & source, std::vector<double> &distance, const int & target = -1);
	// compute path from target to source
	void computePath(const AdjacencyList & graph, const int & target, std::vector<int> & path);

private:

	struct Greater {
		bool operator() (const GraphPair & pair1, const GraphPair & pair2) {
			return pair1.second > pair2.second;
		}
	};
	struct Less {
		bool operator()(const GraphPair & pair1, const GraphPair & pair2) {
			return pair1.second < pair2.second;
		}
	};

	AdjacencyList _distanceGraph; // GraphPair means <index of the node, shortest distance of the node>
};

#endif