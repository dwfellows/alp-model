% Third-order quadrature function
function [int] = quadrature(val1, val2, val3)

  % Weights obtained from Zienkiewics and modified to ensure that transformed surface integral provides A = 1/2.
  w1 = 1/6;
  w2 = 1/6;
  w3 = 1/6;

  int = w1*val1 + w2*val2 + w3*val3;

end

