function blur_sigma = calSigmaValue(Y, Sxf, Syf, dim, tzSquare)

	%[h, w] = size(Y);
	blur_sigma = 1;
	TVperNN = 50000;
	while (TVperNN > 1870)%0.03125 1302.083 1736.2
		blur_sigma = blur_sigma + 1;
		blur_size = 4 * blur_sigma + 1;
		params.blur_kernel  = fspecial('gaussian', blur_size, blur_sigma);
		YBlur = imfilter(Y,params.blur_kernel,'symmetric');
		
		localMax = ordfilt2(YBlur, 9, true(3));
		localMin = ordfilt2(YBlur, 1, true(3));
		Variance = max(localMax - YBlur, YBlur - localMin);
		%TVperNN = mean2(Variance) * focal_length * focal_length / tzSquare;
		TVperNN = mean2(Variance) * (2*Sxf*dim.marker_w) * (2*Syf*dim.marker_h) / tzSquare;
	end
	
end