#ifndef DIJKSTRA_H
#define DIJKSTRA_H

/*
	compute shortest distance and path of a graph represented by adjacency list via Dijkstra algorithm

	todo: 改成1次计算所有距离的，

	jjcao @ 2016
*/
#include <iostream>
#include <queue>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <queue>


typedef std::pair<int, double> GraphPair;
typedef std::vector<std::vector<GraphPair> > AdjacencyList;

class Dijkstra {
public:
    // variables
    AdjacencyList graph;

    std::vector<double> distance;
    std::vector<int> parent;
    std::vector<int> path;
    
    unsigned counter;

	void buildGraph(AdjacencyList &fromMesh);
    
    void initialization(const int & source);
    
    void computeDistances(const int & target);
    
    bool createPath(const int & source, const int & target);
    
    void printPath(const int & target);
    
    // comperator
    struct Comperator {
        int operator() (const GraphPair & pair1, const GraphPair & pair2) {
            return pair1.second > pair2.second;
        }
    };
    
    std::priority_queue<GraphPair, std::vector<GraphPair>, Comperator> priorityQueue;

    
//public:
//	~Dijkstra()
//	{
//		graph.clear();
//		distance.clear();
//		parent.clear();
//		path.clear();
//		counter = 0;
//	}
};

#endif