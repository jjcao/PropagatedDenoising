#include "dijkstra.h"

// license
//void Dijkstra::license() {
//    std::cout << "The MIT License (MIT)\n"
//              << "Copyright (c) 2014 Jonathan Baumann\n"
//              << "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n\n"
//              << "The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\n"
//              <<  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n\n"
//    << std::endl;
//}

void Dijkstra::buildGraph(AdjacencyList &fromMesh)
{
	counter = fromMesh.size();
	graph.resize(counter);
	int i = 0;
	for (std::vector<std::vector<GraphPair> >::iterator iter = fromMesh.begin(); iter != fromMesh.end(); iter++)
	{
		graph[i] = *iter;
		//G_vertices.insert(std::unordered_map<int, const std::unordered_map<int, double>>::value_type((*f_it).idx(), G_edges));
		//for (std::vector<graphPair>::iterator iter_p = (*iter).begin(); iter_p != (*iter).end(); iter_p++)
		//{
		//	//fromMesh[f_it->idx()].push_back(std::vector<graphPair>::value_type(ff_it->idx(), NormalDistance(normals[(*f_it).idx()], normals[(*ff_it).idx()])));
		//	graph[i].push_back(std::make_pair((*iter_p).first, (*iter_p).second));
		i++;
		//}
	}
}
// initalizes the distances and the parent elements
void Dijkstra::initialization(const int & source) {
    
    distance.resize(counter);
    parent.resize(counter);
    
    for(size_t i = 0; i < counter; ++i) {
        //distance[i] = std::numeric_limits<float>::max();
		distance[i] = std::numeric_limits<double>::max();
        parent[i] = -1;
    }
    
    distance[source] = 0.0;
    
    priorityQueue.push(std::make_pair(source, distance[source]));
}

// computes the distances for the nodes
void Dijkstra::computeDistances(const int & target) {
    
    while(!priorityQueue.empty()) {
        int u = priorityQueue.top().first;
        
        if(u == target) {
            break;
        }
        
        priorityQueue.pop();
        
        for(unsigned i = 0; i < graph[u].size(); ++i) {
            int v = graph[u][i].first;
            double w = graph[u][i].second;
            
            if(distance[v] > distance[u] + w) {
                distance[v] = distance[u] + w;
                parent[v] = u;
                
                priorityQueue.push(std::make_pair(v, distance[v]));
            }
        }
    }
}

// searches the best path
bool Dijkstra::createPath(const int & source, const int & target) {
    
    bool hasPath = true;
    
    int p = target;
    
    path.push_back(target);
    
    while(p != source) {
        p = parent[p];
        
        if(p == -1) {
            hasPath = false;
            break;
        }
        
        path.push_back(p);
    }
    
    return hasPath;
}

// prints the path to the console
void Dijkstra::printPath(const int & target) {
    
    std::cout << distance[target] << std::endl;
    
    for(int i = path.size() - 1; i >= 0; --i) {
        std::cout << path[i] << " ";
    }
    
    std::cout << std::endl;
}