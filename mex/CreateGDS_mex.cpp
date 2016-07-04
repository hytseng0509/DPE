
#include "mex.h"
#include <math.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{

	/* Retrieve the input data */
	double *tx_steps = mxGetPr(prhs[0]);
	double *ty_steps = mxGetPr(prhs[1]);
	double *tz_steps = mxGetPr(prhs[2]);
	double *rx_steps = mxGetPr(prhs[3]);
	double *rz0_steps = mxGetPr(prhs[4]);
	double *rz1_steps = mxGetPr(prhs[5]);

	/* Create an mxArray for the output data */
	plhs[0] = mxCreateDoubleMatrix(6, 13, mxREAL );  // 6 for variable for pose variable

	/* Create a pointer to the output data */
	double *poses = mxGetPr(plhs[0]);
	
	// MAIN LOOP
	int gridInd = 0;
  int ind;
  #pragma omp parallel for private(ind) num_threads(8)
	for (ind = 0; ind < 13; ind++) {
		poses[6*ind] = tx_steps[ind];
		poses[6*ind+1] = ty_steps[ind];
		poses[6*ind+2] = tz_steps[ind];
		poses[6*ind+3] = rx_steps[ind];
		poses[6*ind+4] = rz0_steps[ind];
		poses[6*ind+5] = rz1_steps[ind];
	}
}