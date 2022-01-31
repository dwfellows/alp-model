% Triangular interpolation function for data
function [quant] = tri_interp(interp_loc, data);

  quant = data(1)*interp_loc(1) + data(2)*interp_loc(2) + data(3)*interp_loc(3);

end

