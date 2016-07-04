function [ex_mat bestEa] = GDS(marker, img, in_mat, inipose, minDim, delta, bounds, steps, dim, photometricInvariance, verbose)
  
  % smooth the img
  blur_sigma = 2;
  blur_size = 4 * blur_sigma;
  params.blur_kernel  = fspecial('gaussian', blur_size, blur_sigma);
  marker = imfilter(marker, params.blur_kernel, 'symmetric');
  img = imfilter(img, params.blur_kernel, 'symmetric');
  
  % rgb to ycbcr
	marker = rgb2ycbcrNorm(marker);
	img = rgb2ycbcrNorm(img);
  
  % sampling
  epsilon = 0.075;
  numPoints = round(20/epsilon^2);
  xs = randi(dim.marker.w, [1,numPoints]);
  ys = randi(dim.marker.h, [1,numPoints]);

  %% main loop
  deltaFact = 1.511;%1.308;
  level = 1;
  newDelta = delta;
  bestEas = ones(1,3);
  while (1)
    if (verbose)
			fprintf('-- level %d -- delta %f, ', level, newDelta);
		end
    % create 13-points searching pattern
    poses = CreateGDSpattern(inipose, steps, bounds, dim.marker_w, dim.marker_h);

    % pose 2 transformation matrix
    [trans_mex, insiders] = Poses2Trans_mex(poses', int32(dim.img.h), int32(dim.img.w),...
													in_mat(1,1), in_mat(1,3), in_mat(2,2), in_mat(2,3),...
													dim.marker_w, dim.marker_h);
    inBoundaryInds = find(insiders);
    trans_mex = trans_mex(:,inBoundaryInds);
    poses = poses(inBoundaryInds,:);
    if (numel(poses) == 0)
      [~, ex_mat] = getTransAndExMatrix(inipose, in_mat);
      bestEa = 1;
      break;
    end
    
    % evaluate Ea of all poses
		if (photometricInvariance)
			distances = EvaluateEa_photo_mex(marker(:,:,1)',img(:,:,1)',...
											marker(:,:,2)',img(:,:,2)',...
											marker(:,:,3)',img(:,:,3)',trans_mex,int32(xs),int32(ys),dim.marker_w,dim.marker_h);
		else
			distances = EvaluateEa_color_mex(marker(:,:,1)',img(:,:,1)',...
											marker(:,:,2)',img(:,:,2)',...
											marker(:,:,3)',img(:,:,3)',trans_mex,int32(xs),int32(ys),dim.marker_w,dim.marker_h);
		end

    % get best pose with minimun Ea
		[bestEa,ind] = min(distances);
    bestPose = poses(ind,:);
		[~, ex_mat] = getTransAndExMatrix(bestPose, in_mat);    
		if (verbose)
			fprintf('Best Ea %f\n', bestEa);   
		end
		
    % to next level?
		if ((sum(bestPose == inipose) == 6) || (mean(bestEas)*0.999 < bestEa))
      % terminate?
      if (newDelta < 0.005)
        break;
      end
      level = level + 1;
			newDelta = newDelta/deltaFact;
			steps.tx = steps.tx /deltaFact;
			steps.ty = steps.ty /deltaFact;
			steps.tz = steps.tz /deltaFact;
			steps.rx = steps.rx /deltaFact;
			steps.ry = steps.rz0/deltaFact;
			steps.rz = steps.rz1/deltaFact;
			xs = randi(dim.marker.w, [1,numPoints]);
      ys = randi(dim.marker.h, [1,numPoints]);
      bestEas = ones(1,3);
		else
			inipose = bestPose;
      bestEas = [bestEas(2:3) bestEa];
		end
  end
end

function img = rgb2ycbcrNorm(img)
	img = rgb2ycbcr(img);
	img(:,:,1) = (img(:,:,1) - 16/255) / (235/255-16/255);
	img(:,:,2) = (img(:,:,2) - 16/255) / (240/255-16/255);
	img(:,:,3) = (img(:,:,3) - 16/255) / (240/255-16/255);
end