function [corner_x, corner_y] = draw_coordinate(extrinsic, intrinsic, minDim)

  M_4_4 = intrinsic * [extrinsic; 0 0 0 1];     %%    transform marker coor to camera screen
	boundary_x = minDim .* [0 1 0 0]; 
	boundary_y = minDim .* [0 0 1 0];
	boundary_z = minDim .* [0 0 0 1];
	for i = 1:4
		corner = M_4_4*[boundary_x(i); boundary_y(i); boundary_z(i); 1];
		corner_x(i) = corner(1)/corner(3);
		corner_y(i) = corner(2)/corner(3);
	end
end