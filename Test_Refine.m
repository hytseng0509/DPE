function [ex_mat ex_mats] = Test_Refine(ex_mat_ini, M, I, in_mat, minDim, minTz, maxTz, delta, photometricInvariance, needcompile, verbose)
% Test_Refine: test the proposed refinement scheme
%
% Input:
%     - ex_mat_ini: initial extrinsic matrix 3*4
%     - M: target image (double)
%     - I: camera image (double)
%     - in_mat: camera intrinsic matrix 4*4
%     - minDim: length of the shorter side of the target
%     - minTz:  minimum distance between camera and target
%     - maxTz:  maximum distance between camera and target
%     - delta: initial precision parameter (default 0.25)
%     - photometricInvariance: need to be photometric invariant
%     - needcompile: need to compile the mex file
%     - verbose: show the state of the method
% Output:
%     - ex_mat: estimated extrinsic matrix
  AddPaths;
  if (~exist('needcompile','var'))
    needcompile = 0;
  end
  if (needcompile == 1)
    CompileMex;
  end
  
  % preCalculation
  [marker, img, bounds, steps, dim] = preCal(in_mat, M, I, minDim, minTz, maxTz, delta, verbose);
  
  % refine
  [ex_mat ex_mats] = Refine(M, I, in_mat, ex_mat_ini, minDim, delta, bounds, steps, dim, photometricInvariance, verbose);
  
end

