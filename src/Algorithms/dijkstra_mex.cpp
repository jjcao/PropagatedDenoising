/*=================================================================
*
*
* usage:
D = dijkstra_mex(W, startIdx);
* 'W' is the weight matrix. W(i,j) is the cost of moving from i to j,
* now: W is not a sparse matrix.
*
* startIdx is the index of the start point
* D(j) will be distance between the start point and j.
*
* JJCAO, 2016
*
*=================================================================*/
#include <mex.h>
#include "dijkstra.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	// input
	int row = mxGetM(prhs[0]);
	int col = mxGetN(prhs[0]);
	if (col != row)
		mexErrMsgTxt("W should be a square matrix!");

	double *w = mxGetPr(prhs[0]);
	int* startIdx = (int*)mxGetPr(prhs[1]);

	// output
	plhs[0] = mxCreateDoubleMatrix(row, 1, mxREAL);
	double *d = mxGetPr(plhs[0]);

	// compute
	Dijkstra dij;
	AdjacencyList localGraph;
	localGraph.resize(row);

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			if(w[i+j*row]>0)
				localGraph[i].push_back(std::make_pair(j, w[i + j*row]));
		}		
	}

	std::vector<double> distance;
	dij.computeDistances(localGraph, *startIdx, distance);
	for (int i = 0; i < row; ++i)
	{
		d[i] = distance[i];
	}
}