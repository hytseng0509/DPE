function ex_mat = Test_DPE(Marker, img, in_mat, minDim, minTz, maxTz, delta, photometricInvariance, needcompile, verbose)
% Test_DPE: test the proposed 3D pose estimation method
%
% Input:
%     - Marker: target image (double)
%     - img:    camera image (double)
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
  
  [ex_mat, ex_mats] = DPE(Marker,img,in_mat,minDim,photometricInvariance,minTz,maxTz,delta,verbose);
end

