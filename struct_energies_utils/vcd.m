function [dLidxj] = vcd(i,j,coeffs)
  % vcd -- Volume coordinate derivative
  % L_i = (1/6V)(a_i + b_i x + c_i y + d_i z)
  % dL_i / dx_j = (1/6V) e_j; e_1 = b_i, e_2 = c_i, e_3 = d_i

  V = coeffs(1);
  a1 = coeffs(2);
  b1 = coeffs(3);
  c1 = coeffs(4);
  d1 = coeffs(5);
  a2 = coeffs(6);
  b2 = coeffs(7);
  c2 = coeffs(8);
  d2 = coeffs(9);
  a3 = coeffs(10);
  b3 = coeffs(11);
  c3 = coeffs(12);
  d3 = coeffs(13);
  a4 = coeffs(14);
  b4 = coeffs(15);
  c4 = coeffs(16);
  d4 = coeffs(17);

  if (i==1)  %Return dL1 / dx_j
    if (j==1)
      dLidxj = b1;
    elseif (j==2)
      dLidxj = c1;
    elseif (j==3)
      dLidxj = d1;
    end
  elseif (i==2) % Return dL2 / dx_j
    if (j==1)
      dLidxj = b2;
    elseif (j==2)
      dLidxj = c2;
    elseif (j==3)
      dLidxj = d2;
    end
  elseif (i==3) % Return dL3 / dx_j
    if (j==1)
      dLidxj = b3;
    elseif (j==2)
      dLidxj = c3;
    elseif (j==3)
      dLidxj = d3;
    end
  else % Return dL4 / dx_j
    if (j==1)
      dLidxj = b4;
    elseif (j==2)
      dLidxj = c4;
    elseif (j==3)
      dLidxj = d4;
    end
  end

  % Augment by (1/(6V))
  dLidxj = dLidxj / (6*V);

end

