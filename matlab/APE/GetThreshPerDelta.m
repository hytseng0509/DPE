function thresh = GetThreshPerDelta(delta)
	safety = 0.002;
	%if delta > 0.02	
	%	thresh = 0.1803 * delta + 0.0219 - safety;
	%else
	%	thresh = 0.8373 * delta;% + 0.0018 - safety;
	%end
	
	%thresh = 0.2069 * delta + 0.0174 - safety;
	thresh = 0.1869 * delta + 0.0161 - safety;
end