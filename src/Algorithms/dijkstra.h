#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <iostream>
#include <queue>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <queue>

typedef std::pair<int, double> graphPair;
typedef std::vector<std::vector<graphPair> > dijkstraGraph;

class Dijkstra {
public:
    // variables
    dijkstraGraph graph;

    std::vector<double> distance;
    std::vector<int> parent;
    std::vector<int> path;
    
    unsigned counter;
    
    // methods
    //void program(int argc, char * argv[]);
    
    //void license();
    
    //void readTGF(const std::string & filename);
	void buildGraph(std::vector<std::vector<graphPair> > &fromMesh);
    
    void initialization(const int & source);
    
    void computeDistances(const int & target);
    
    bool createPath(const int & source, const int & target);
    
    void printPath(const int & target);
    
    // comperator
    struct Comperator {
        int operator() (const graphPair & pair1, const graphPair & pair2) {
            return pair1.second > pair2.second;
        }
    };
    
    std::priority_queue<graphPair, std::vector<graphPair>, Comperator> priorityQueue;

    
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