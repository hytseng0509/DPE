function [ex_mat ex_mats] = Refine(marker, img, in_mat, ex_mat, minDim, delta, bounds, steps, dim, photometricInvariance, verbose)
  
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