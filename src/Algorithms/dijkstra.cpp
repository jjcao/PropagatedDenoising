#include "dijkstra.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>

void Dijkstra::computeDistances(const AdjacencyList & graph, const int & source, std::vector<double> &distance, const int & target)
{
	// initialize
	int numNode(graph.size());
	distance.clear();
	distance.resize(numNode);

	for (size_t i = 0; i < numNode; ++i)
		distance[i] = std::numeric_limits<double>::max();

	distance[source] = 0.0;

	std::priority_queue<GraphPair, std::vector<GraphPair>, Greater> priorityQueue;
	priorityQueue.push(std::make_pair(source, distance[source])); 
	

	// compute    
    while(!priorityQueue.empty()) 
	{
        int u = priorityQueue.top().first;        
        if(u == target)    break;        
        priorityQueue.pop();
        
        for(unsigned i = 0; i < graph[u].size(); ++i) 
		{
            int v = graph[u][i].first;
            double w = graph[u][i].second;
            
            if(distance[u] + w < distance[v])
			{
				distance[v] = distance[u] + w;                
                priorityQueue.push(std::make_pair(v, distance[v]));
            }
        }
    }

	// save it for compute shortest path later.
	_distanceGraph.clear();
	_distanceGraph.resize(numNode);
	for (int i = 0; i < numNode; ++i)
	{
		for (int j = 0; j < graph[i].size(); ++j)
		{
			_distanceGraph[i].push_back(std::make_pair(graph[i][j].first, distance[graph[i][j].first]));
		}		
	}
	
}

void Dijkstra::computePath(const AdjacencyList & graph, const int & target, std::vector<int> & path)
{
	int numNode(graph.size());
	if (target >= numNode) return;

	path.clear();
	path.push_back(target);
	int u = target;
	double distBefore(std::numeric_limits<double>::max());
	while (true)
	{
		std::vector<GraphPair> & tmpList = _distanceGraph[u];
		std::vector<GraphPair>::iterator it = std::min_element(tmpList.begin(), tmpList.end(), Less());
		if (it->second < distBefore)
			distBefore = it->second; 
		else // source is achieved or there is node which is not visited when computeDistance();
			return;			

		u = it->first;
		path.push_back(u);
	}
}