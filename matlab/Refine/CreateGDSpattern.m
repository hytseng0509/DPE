function poses = CreateGDSpattern(inipose, steps, bounds, marker_w, marker_h);
  
  weight = inipose(3) + norm([marker_w, marker_h])*sin(inipose(4));
  tx_steps = zeros(1,13) + inipose(1);
  tx_steps(2) = tx_steps(2) + steps.tx*weight;
  tx_steps(3) = tx_steps(3) - steps.tx*weight;
  
  ty_steps = zeros(1,13) + inipose(2);
  ty_steps(4) = ty_steps(4) + steps.ty*weight;
  ty_steps(5) = ty_steps(5) - steps.ty*weight;
  
  tz_steps = zeros(1,13) + inipose(3);
  tz_steps(6) = min(tz_steps(6) + (tz_steps(6))^2 * steps.tz / (1 - steps.tz*tz_steps(6)), bounds.tz(2));
  tz_steps(7) = max(tz_steps(7) - (tz_steps(7))^2 * steps.tz / (1 + steps.tz*tz_steps(7)), bounds.tz(1));
  
  rx_steps = zeros(1,13) + inipose(4);
  positve  = asin(2 - 1/(1/(2-sin(inipose(4))) + steps.rx));
  negative = asin(2 - 1/(1/(2-sin(inipose(4))) - steps.rx));
  if (imag(positve)~=0)
    positve = bounds.rx(2);
  end
  if (imag(negative)~=0)
    negative = bounds.rx(1);
  end
  rx_steps(8) = min(positve, bounds.rx(2));
  rx_steps(9) = max(negative, bounds.rx(1));

  rz0_steps = zeros(1,13) + inipose(5);
  rz0_steps(10) = rz0_steps(10) + steps.rz0;
  rz0_steps(11) = rz0_steps(11) - steps.rz0;
  
  rz1_steps = zeros(1,13) + inipose(6);
  rz1_steps(12) = rz1_steps(10) + steps.rz1;
  rz1_steps(13) = rz1_steps(11) - steps.rz1;

  poses_mex = CreateGDS_mex(tx_steps, ty_steps, tz_steps, rx_steps, rz0_steps, rz1_steps);
  poses = poses_mex';
  
end