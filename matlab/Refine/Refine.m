function [ex_mat ex_mats] = Refine(marker, img, in_mat, ex_mat, minDim, delta, bounds, steps, dim, photometricInvariance, verbose)
  
  % smooth images
	blur_size = 4 * 2 + 1;
	params.blur_kernel  = fspecial('gaussian', blur_size, 2);
	marker = imfilter(marker,params.blur_kernel,'symmetric');
	img = imfilter(img,params.blur_kernel,'symmetric');

	% rgb to ycbcr
	marker = rgb2ycbcrNorm(marker);
	img = rgb2ycbcrNorm(img);
  
  % get 2 poses
  poses = get2Poses(in_mat, ex_mat, minDim);
  
  % GDS refinement
  ex_mats = zeros(3, 4, 2);
  Eas = ones(1, 2);
  for i = 1: size(poses, 1)
    [ex_mats(:, :, i) Eas(i)] = GDS(marker, img, in_mat, poses(i, :), minDim, delta, bounds, steps, dim, photometricInvariance, verbose);
  end
  
  % get the one with smaller Ea
  if (Eas(1) <= Eas(2))
    ex_mat = ex_mats(:, :, 1);
  else
    if (verbose)
      fprintf('switch to the second pose!\n');
    end
    ex_mat = ex_mats(:, :, 2);
  end
  
end

function img = rgb2ycbcrNorm(img)
	img = rgb2ycbcr(img);
	img(:,:,1) = (img(:,:,1) - 16/255) / (235/255-16/255);
	img(:,:,2) = (img(:,:,2) - 16/255) / (240/255-16/255);
	img(:,:,3) = (img(:,:,3) - 16/255) / (240/255-16/255);
end