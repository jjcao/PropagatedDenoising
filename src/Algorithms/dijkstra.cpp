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

// reads the tgf file and builds the graph
//void Dijkstra::readTGF(const std::string & filename) {
//    
//    std::ifstream tgfFile(filename);
//    
//    if(!tgfFile.is_open()) {
//        std::cerr << "Could not open the file: " << filename << std::endl;
//    }
//    
//    // which part of the file; false = first part; true = second part
//    bool part = false;
//    
//    // linebuffer
//    std::string linebuffer;
//    
//    // count the elements of the file
//    counter = 0;
//    
//    // read every line
//    while(std::getline(tgfFile, linebuffer)) {
//        
//        // if a # is read
//        if(linebuffer.substr(0, 1) == "#") {
//            part = true;
//            graph.resize(counter);
//        }
//        
//        else {
//            std::istringstream split(linebuffer);
//            std::vector<std::string> tokens;
//            
//            if(part) {
//                // save the lines in the graph
//                for(std::string each; std::getline(split, each, ' '); tokens.push_back(each));
//                
//                int id1 = atoi(tokens[0].c_str());
//                int id2 = atoi(tokens[1].c_str());
//                float weight = atof(tokens[2].c_str());
//                
//                graph[id1].push_back(std::make_pair(id2, weight));
//            }
//            
//            else {
//                ++counter;
//            }
//        }
//    }
//    
//    tgfFile.close();
//}
void Dijkstra::buildGraph(std::vector<std::vector<graphPair> > &fromMesh)
{
	counter = fromMesh.size();
	graph.resize(counter);
	int i = 0;
	for (std::vector<std::vector<graphPair> >::iterator iter = fromMesh.begin(); iter != fromMesh.end(); iter++)
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

// runs the program and calls the functions
//void Dijkstra::program(int argc, char * argv[]) {
//    
//    // argc = 4 because the argv at position 0 is the program itself
//    if(argc != 4) {
//        std::cout << "Dijkstra needs 3 parameters:\n     inputfile in tgf format\n     source\n     target";
//    }
//    
//    std::string filename = argv[1];
//    int source = atoi(argv[2]);
//    int target = atoi(argv[3]);
//    
//    // license
//    license();
//    
//    // read the tgf file
//    readTGF(filename);
//    
//    // initalization of distance and parent vector
//    initialization(source);
//    
//    // compute the distances (dijkstra)
//    computeDistances(target);
//    
//    // create the path and prints it out
//    if(createPath(source, target)) {
//        printPath(target);
//    }
//    
//    else {
//        std::cout << "no path" << std::endl;
//    }
//    
//}