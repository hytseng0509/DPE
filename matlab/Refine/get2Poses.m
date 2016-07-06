function poses = get2Poses(in_mat, ex_mat, minDim);
  
  poses = [];
  tgt = [-minDim,minDim,minDim,-minDim;-minDim,-minDim,minDim,minDim;0,0,0,0;1,1,1,1];
  src = ex_mat*tgt;
  src(1,:) = src(1,:)./src(3,:);
  src(2,:) = src(2,:)./src(3,:);
  [RR, tt, error, flag] = OPnP(tgt(1:3,:), src(1:2,:));
  if (flag)
    poses = exMat2Rz0RxRz1(ex_mat(:,1:3), ex_mat(:,4));
  else
    numPoses = min(size(RR, 3), 2);
    poses = zeros(numPoses, 6);
    for i = 1:numPoses
      poses(i, :) = exMat2Rz0RxRz1(RR(:, :, i), tt(:, i));
    end
  end
end

function pose = exMat2Rz0RxRz1(R, t)
	rx = acos(R(3,3)) - pi;
	if (rx == 0)
		rz0 = 0;
		rz1 = atan2(-R(1,2), R(1,1));
	else
		rz0 = atan2(R(1,3), -R(2,3));
		rz1 = atan2(R(3,1), R(3,2));
	end
	pose = [t(1) t(2) t(3) rx rz0 rz1];
end