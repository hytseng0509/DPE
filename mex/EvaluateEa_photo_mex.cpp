#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
  
  // parameters
  int h1, w1, h2, w2;
  int r1x, r1y, r2x, r2y;
  int numPoints, numPoses;
  
  // input variables
  double *img1, *img1_cr, *img1_cb, *img2, *img2_cr, *img2_cb;
  double *trans;     // element = 16
  int *xs, *ys;
  
  // helper variables
  double *xs_centered, *ys_centered;
  double *valsI1, *valsI1_cr, *valsI1_cb;

  // output variables
  double *distances;
  
  /* Find the dimensions of the data */
  h1 = mxGetN(prhs[0]);
  w1 = mxGetM(prhs[0]);
  h2 = mxGetN(prhs[1]);
  w2 = mxGetM(prhs[1]);
  numPoses = mxGetN(prhs[6]);
  numPoints = mxGetN(prhs[7]);
  r1x = 0.5*(w1-1);
  r1y = 0.5*(h1-1);
  r2x = 0.5*(w2-1);
  r2y = 0.5*(h2-1);

  /* Create an mxArray for the output data */
  plhs[0] = mxCreateDoubleMatrix(1, numPoses, mxREAL );
  
  /* Create an mxArrays for temporary data */
  xs_centered = (double *)malloc(numPoints*sizeof(double));
  ys_centered = (double *)malloc(numPoints*sizeof(double));
  valsI1 = (double *)malloc(numPoints*sizeof(double));
	valsI1_cr = (double *)malloc(numPoints*sizeof(double));
	valsI1_cb = (double *)malloc(numPoints*sizeof(double));
   
  /* Retrieve the input data */
  img1 = mxGetPr(prhs[0]);
  double* tmp_img2 = mxGetPr(prhs[1]);
	img1_cb = mxGetPr(prhs[2]);
	double* tmp_img2_cb = mxGetPr(prhs[3]);
	img1_cr = mxGetPr(prhs[4]);
	double* tmp_img2_cr = mxGetPr(prhs[5]);
  trans = mxGetPr(prhs[6]);
  xs = (int*)mxGetPr(prhs[7]);
  ys = (int*)mxGetPr(prhs[8]);
	double marker_w = mxGetScalar(prhs[9]);
	double marker_h = mxGetScalar(prhs[10]);
    
  //img2 is of height 3*img2 (this padding is for not needing to check bounds)
  img2 = (double*)malloc(5*h2*w2*sizeof(double));
  memset(img2, 2, 5*h2*w2*sizeof(double));
  memcpy(img2+2*h2*w2, tmp_img2, h2*w2*sizeof(double));
  
  img2_cr = (double*)malloc(5*h2*w2*sizeof(double));
  memset(img2_cr, 2, 5*h2*w2*sizeof(double));
  memcpy(img2_cr+2*h2*w2, tmp_img2_cr, h2*w2*sizeof(double));
	
	img2_cb = (double*)malloc(5*h2*w2*sizeof(double));
  memset(img2_cb, 2, 5*h2*w2*sizeof(double));
  memcpy(img2_cb+2*h2*w2, tmp_img2_cb, h2*w2*sizeof(double));
    
  /*Centered pointes for marker [-0.5 0.5]*/
  for (int i = 0 ; i < numPoints ; i++) {
    xs_centered[i] =  double(xs[i]-(r1x+1))/(w1*0.5)*marker_w;
    ys_centered[i] = -double(ys[i]-(r1y+1))/(h1*0.5)*marker_h;
  }
    
  /*Precalculating source point indices into I1 (and the values themselves)*/
  for (int j = 0; j < numPoints ; j++) {
    valsI1[j] = img1[(ys[j] - 1)*w1 + xs[j]-1]; // -1 is for c
    valsI1_cr[j] = img1_cr[(ys[j] - 1)*w1 + xs[j]-1]; // -1 is for c
    valsI1_cb[j] = img1_cb[(ys[j] - 1)*w1 + xs[j]-1]; // -1 is for c
  }
    
    
  /* Create a pointer to the output data */
  distances = mxGetPr(plhs[0]);


      
  // MAIN LOOP
  int i;
  #pragma omp parallel for private(i) num_threads(8)
  for (i = 0 ; i < numPoses ; i++) {       
    int targetPoint_x, targetPoint_y;
    int targetInd;
    double a11, a12, a13, a21, a22, a23, c1, c2, c3;
    a11 = trans[16*i];
    a12 = trans[16*i+1];
    a13 = trans[16*i+3];
    a21 = trans[16*i+4];
    a22 = trans[16*i+5];
    a23 = trans[16*i+7];
    c1  = trans[16*i+8];
    c2  = trans[16*i+9];
    c3  = trans[16*i+11];
    
    double score = 0;
    
    double* ptrVals = valsI1;
    double* ptrXsc = xs_centered;
    double* ptrYsc = ys_centered;
    
    double sumXi =0; double sumYi =0;
    double sumXiSqrd =0; double sumYiSqrd =0; double Xi, Yi;

    double* ys_target = (double *)malloc(numPoints*sizeof(double));
    double* ys_target_cr = (double *)malloc(numPoints*sizeof(double));
    double* ys_target_cb = (double *)malloc(numPoints*sizeof(double));
    for (int j = 0; j < numPoints ; j++) {
			double marker_x = double(*ptrXsc);
			double marker_y = double(*ptrYsc);
      targetPoint_x = int(double(a11*marker_x  + a12*marker_y + a13) / double(c1*marker_x  + c2*marker_y + c3) + 0.5); // includes rounding(perspective form)
      targetPoint_y = int(double(a21*marker_x  + a22*marker_y + a23) / double(c1*marker_x  + c2*marker_y + c3) + 0.5 + 2*h2); // includes rounding
      targetInd = (targetPoint_y - 1)*w2 + targetPoint_x - 1; // -1 is for c
      
      Xi = (*ptrVals) ;
      Yi = img2[targetInd];
      ys_target[j] = Yi;
			ys_target_cr[j] = img2_cr[targetInd];
			ys_target_cb[j] = img2_cb[targetInd];
      
      sumXi += Xi;
      sumYi += Yi;
      sumXiSqrd += (Xi*Xi);
      sumYiSqrd += (Yi*Yi);
      ptrVals++;
      ptrXsc++;
      ptrYsc++;
    }

    double epsilon = 0.0000001;
    double sigX = sqrt((sumXiSqrd-(sumXi*sumXi)/numPoints)/numPoints) + epsilon;
    double sigY = sqrt((sumYiSqrd-(sumYi*sumYi)/numPoints)/numPoints) + epsilon;
    double meanX = sumXi/numPoints;
    double meanY = sumYi/numPoints;
    double sigXoversigY = sigX/sigY;
          
    // Variable that stores a sum used repeatadly in the computation: -meanX+sigXoversigY*meanY
    double faster = -meanX+sigXoversigY*meanY;
    for (int j = 0; j < numPoints ; j++) {
      score += (0.5*fabs(valsI1[j]-sigXoversigY*ys_target[j]+faster)+0.25*fabs(valsI1_cr[j]-ys_target_cr[j])+0.25*fabs(valsI1_cb[j]-ys_target_cb[j]));
    }

    free(ys_target);
    free(ys_target_cr);
    free(ys_target_cb);
      
    distances[i] = score/numPoints;
  }

  /* Free the allocated arrays */
  free(xs_centered);
  free(ys_centered);
  free(valsI1);
	free(valsI1_cr);
	free(valsI1_cb);
  free(img2);
	free(img2_cb);
	free(img2_cr);
}