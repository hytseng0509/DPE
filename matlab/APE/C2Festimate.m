function [ex_mat,newDelta,steps] = C2Festimate(marker, img, in_mat, bounds, steps, dim, epsilon, delta, photometricInvariance, verbose)
	
	% get random sample of pixels for computing Ea
	numPoints = round(10/epsilon^2);
	xs = randi(dim.marker.w,[1,numPoints]);
	ys = randi(dim.marker.h,[1,numPoints]);
	
	% create the epsilon-cover pose set
	[poses, numPoses] = CreateEpsilonCoverSet(in_mat, bounds, steps, dim);
	if (numPoses > 71000000)
		error('more than 71 million poses!');
	end
	
	% start C2F estimate
	deltaFact = 1.511;
	level = 0;
	bestEa = [];
  bestDists = zeros(1 ,8);
	newDelta = delta;
	totTime = 0;
  if (photometricInvariance)
    c1 = 0.075; c2 = 0.15;
  else
    c1 = 0.05; c2 = 0.1;
  end
	while (1)
		level = level + 1;
		if (verbose)
			fprintf('-- level %d -- delta %f, ', level, newDelta);
		end
    
		% calculate transformation matrix from poses
		Poses2TransST = tic;
		[trans_mex, insiders] = Poses2Trans_mex(poses', int32(dim.img.h), int32(dim.img.w),...
												in_mat(1,1), in_mat(1,3), in_mat(2,2), in_mat(2,3),...
												dim.marker_w, dim.marker_h);
		inBoundaryInds = find(insiders);
		trans_mex = trans_mex(:,inBoundaryInds);
		poses = poses(inBoundaryInds,:);
    numPoses = size(poses, 1);
		Poses2TransTime = toc(Poses2TransST);
    if (verbose)
			fprintf('Number of Poses %d, ', numPoses);
		end
    
		% evaluate Ea of all poses
		EvaluateEaST = tic;
		if (photometricInvariance)
			distances = EvaluateEa_photo_mex(marker(:,:,1)',img(:,:,1)',...
											marker(:,:,2)',img(:,:,2)',...
											marker(:,:,3)',img(:,:,3)',trans_mex,int32(xs),int32(ys),dim.marker_w,dim.marker_h);
		else
			distances = EvaluateEa_color_mex(marker(:,:,1)',img(:,:,1)',...
											marker(:,:,2)',img(:,:,2)',...
											marker(:,:,3)',img(:,:,3)',trans_mex,int32(xs),int32(ys),dim.marker_w,dim.marker_h);
		end
		EvaluateEaTime = toc(EvaluateEaST);
		totTime = totTime + Poses2TransTime + EvaluateEaTime;
		
		% get best pose with minimun Ea
		[bestEa,ind] = min(distances);
		bestConfig = poses(ind,:);
		[~, ex_mat] = getTransAndExMatrix(poses(ind,:), in_mat);    
		bestDists(level) = bestEa;
		if (verbose)
			fprintf('Best Ea %f, ', bestEa);   
		end
		
		% terminate or not
		if ( (bestEa < 0.005) || ((level > 4) && (bestEa < 0.015)) || ...
				((level > 3) && (bestEa > mean(bestDists(level-3:level-1))*0.97)) || ...
				(level > 7) )
			break
		end
		
		% select poses within threshold to be in the next round
		[~, goodPoses, percentage, tooHighPercentage] = GetPosesByDistance(poses,bestEa,newDelta,distances,verbose);
		if (verbose)
			fprintf('Survived percentage %f\n',percentage * 100);
		end
		
		% expand the pose set for next round
		%if (photometricInvariance == 1)
		%	constraint1 = 10^7;
		%	constraint2 = 7.5*10^6;
		%else
		%	constraint1 = 7.5*10^6;
		%	constraint2 = 5*10^6;
		%end
    %
		%if ((tooHighPercentage && (bestEa > 0.05) && ((level==1) && (size(poses, 1) < constraint1)) ) || ...  %%7.5
		%	(                     (bestEa > 0.10) && ((level==1) && (size(poses, 1) < constraint2)) ) )  %% 5*10^6
		if ((level==1) && ((tooHighPercentage && (bestEa > c1) && (numPoses < 7500000)) || ((bestEa > c2) && (numPoses < 5000000)) ) )
      fact = 0.9;
			newDelta = newDelta*fact;
			level = 0;
			steps.tx = fact*steps.tx;
			steps.ty = fact*steps.ty;
			steps.tz = fact*steps.tz;
			steps.rx = fact*steps.rx;
			steps.rz0 = fact*steps.rz0;
			steps.rz1 = fact*steps.rz1;
			[poses, ~] = CreateEpsilonCoverSet(in_mat, bounds, steps, dim);
		else
			prevDelta = newDelta;
			newDelta = newDelta/deltaFact;
			[expandedPoses, steps] = ExpandPoses(goodPoses,steps,80,deltaFact,bounds,dim.marker_w,dim.marker_h);
			poses = [goodPoses ; expandedPoses];
		end
		
		% re-sample pixels to calculate Ea
		xs = randi(dim.marker.w,[1,numPoints]);
		ys = randi(dim.marker.h,[1,numPoints]);
	end
	
	
	% for output
	sampledError = bestEa;
	return;
end

function [thresh,goodPoses,percentage,tooHighPercentage] = GetPosesByDistance(poses,bestEa,newDelta,distances,verbose)

	thresh = bestEa + GetThreshPerDelta(newDelta);
	goodPoses = poses(distances <= thresh, :);
	numPoses = size(goodPoses,1);
	percentage = numPoses/size(poses,1);
	while (numPoses > 27000)
		thresh = thresh * 0.99;
		goodPoses = poses(distances <= thresh, :);
		numPoses = size(goodPoses,1);
		%if (verbose)
		%	fprintf('too many survivals...\n');
		%end
	end
	
	tooHighPercentage = (percentage > 0.1);

end





