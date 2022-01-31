function [fluid_node_loc_barycentric] = cart2bary(xn, xf)
  % Return location of fluid node in barycentric (volume) coordinates
  % xn: 10x3 array of 10-node tetrahedron coordinates in Cartesian space
  % xf: 1x3 array of fluid node location in Cartesian space

  xi0 = [(1/4)*rand(1); (1/4)*rand(1); (1/4)*rand(1)];
%  xi0 = [0.25; 0.25; 0.25];
  options = optimset('tolFun', 1e-12, 'display', 'iter');
  %options = optimset('tolFun', 1e-12, 'display', 'off');
  xi_solve = fsolve(@(xi) nonlinearBarySystem(xi,xn,xf),xi0,options);

  fluid_node_loc_barycentric = [1 - xi_solve(1) - xi_solve(2) - xi_solve(3); xi_solve(1); xi_solve(2); xi_solve(3)];

end
