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
  double tz_min  = mxGetScalar(prhs[0]);
  double tz_max  = mxGetScalar(prhs[1]);
  double rx_min  = mxGetScalar(prhs[2]);
  double rx_max  = mxGetScalar(prhs[3]);
  double rz_min  = mxGetScalar(prhs[4]);
  double rz_max  = mxGetScalar(prhs[5]);

	double tx_step = mxGetScalar(prhs[6]);
  double ty_step = mxGetScalar(prhs[7]);
	double tz_step = mxGetScalar(prhs[8]);
	double rx_step = mxGetScalar(prhs[9]);
	double rz0_step = mxGetScalar(prhs[10]);
	double rz1_step = mxGetScalar(prhs[11]);
  double Sxf = mxGetScalar(prhs[12]);
  double img_w = mxGetScalar(prhs[13]);
  double Syf = mxGetScalar(prhs[14]);
  double img_h = mxGetScalar(prhs[15]);
	double marker_w = mxGetScalar(prhs[16]);
	double marker_h = mxGetScalar(prhs[17]);

  double tz_mid = sqrt(tz_max*tz_min);
    
  // Pre calculate size
	double tz = tz_min;
  int numTranslate = 0;
	double length = sqrt(marker_w*marker_w + marker_h*marker_h);
  while (tz <= tz_max) {
    double rz1 = rz_min;
		while (rz1 <= rz_max) {
      double rx = rx_min;
			while (rx >= -rx_max) {
        double rz0 = rz_min;
				while (rz0 <= rz_max) {
          double tx_w = img_w * tz / Sxf - marker_h;
          double tx = -tx_w;
					while (tx <= tx_w) {
            double ty_w = img_h * tz / Syf - marker_h;
            double ty = -ty_w;
						while (ty < ty_w){
							numTranslate++;
              ty += ty_step * (tz+length*sin(rx));
						}
						tx += tx_step * (tz+length*sin(rx));
					}
					rz0 += rz0_step;
          if (rx == 0)
            rz0 = rz_max+1;
				}
        double sinValuey = 1 / (1/(2+sin(rx)) + rx_step) - 2;
        if (sinValuey <= 1 && sinValuey >= -1)
          rx = asin(sinValuey);
        else
          rx = -rx_max - 1;
			}
      rz1 += rz1_step;
		}
    tz += tz*tz * tz_step / (1 - tz_step*tz);
	}

  plhs[0] = mxCreateDoubleMatrix(6, numTranslate, mxREAL );  // 6 for variable for pose variable
  double *configs = mxGetPr(plhs[0]);
    
	// MAIN LOOP
	int gridInd = 0;
  tz = tz_min;
  while (tz <= tz_max) {
    double rz1 = rz_min;
		while (rz1 <= rz_max) {
      double rx = rx_min;
			while (rx >= -rx_max) {
        double rz0 = rz_min;
				while (rz0 <= rz_max) {
          double tx_w = img_w * tz / Sxf - marker_h;
          double tx = -tx_w;
					while (tx <= tx_w) {
            double ty_w = img_h * tz / Syf - marker_h;
            double ty = -ty_w;
						while (ty < ty_w) {
              configs[gridInd++] = tx;
              configs[gridInd++] = ty;
              configs[gridInd++] = tz;
              configs[gridInd++] = -rx;
              configs[gridInd++] = rz0;
              configs[gridInd++] = rz1;
              ty += ty_step * (tz+length*sin(rx));
            }
            tx += tx_step * (tz+length*sin(rx));
          }
          rz0 += rz0_step;
          if (rx == 0)
            rz0 = rz_max+1;
				}
        double sinValuey = 1 / (1/(2+sin(rx)) + rx_step) - 2;
        if (sinValuey <= 1 && sinValuey >= -1)
          rx = asin(sinValuey);
        else
          rx = -rx_max - 1;
			}
      rz1 += rz1_step;
		}
    tz += tz*tz * tz_step / (1 - tz_step*tz);
	}
  
}
