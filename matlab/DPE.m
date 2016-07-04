function [ex_mat ex_mats] = DPE(M, I, in_mat, minDim,... % mandatory
   photometricInvariance,minTz,maxTz,delta,verbose,epsilon) %optional                                                                        
	
	% set default values for optional variables
	if (~exist('photometricInvariance','var'))
		photometricInvariance = 0;
	end
	if (~exist('minTz','var'))
		minTz = 3;
	end
	if (~exist('maxTz','var'))
		maxTz = 8;
	end
  if (~exist('delta','var'))
		delta = 0.25;
	end
	if (~exist('verbose','var'))
		verbose = 0;
	end
	if (~exist('epsilon','var'))
		epsilon = 0.15;
	end
	if (~exist('in_mat', 'var'))  %% get intrinsic based on camera image
		[h2,w2,~] = size(img);
		focal_length = norm([h2, w2]);
		in_mat = [focal_length,0,h2/2+0.5,0;0,focal_length,w2/2+0.5,0;0,0,1,0;0,0,0,1];
	end
	
	% ensure the data type of images is double
	if ( ~strcmp(class(I),'double') || ~strcmp(class(M),'double')) %#ok<STISA>
		error('img and marker should both be of class ''double'' (in the range [0,1])');
	end
	
	% preCalculation
  t1 = tic;
	[marker, img, bounds, steps, dim] = preCal(in_mat, M, I, minDim, minTz, maxTz, delta, verbose);
  if (verbose)
    fprintf('pre-time: %f\n', toc(t1));
  end
  
	% coarse-to-fine estimation
  t2 = tic;
	[ex_mat,delta,steps] = C2Festimate(marker,img,in_mat,bounds,steps,dim,epsilon,delta,photometricInvariance,verbose);
	if (verbose)
    fprintf('\npost-time: %f\n', toc(t2));
  end
  
	% refine
  t3 = tic;
  [ex_mat ex_mats] = Refine(M, I, in_mat, ex_mat, minDim, delta, bounds, steps, dim, photometricInvariance, verbose);
  if (verbose)
    fprintf('refine-time: %f\n', toc(t3));
  end
end
