function [interpolation_val] = shape_interpolation(tet_loc, nodal_displacement)

  % Pull off tet_locs
  L1 = tet_loc(1);
  L2 = tet_loc(2);
  L3 = tet_loc(3);
  L4 = tet_loc(4);

  % Pull off nodal displacements
  uI = nodal_displacement(1,:);
  uJ = nodal_displacement(2,:);
  uK = nodal_displacement(3,:);
  uL = nodal_displacement(4,:);
  uM = nodal_displacement(5,:);
  uN = nodal_displacement(6,:);
  uO = nodal_displacement(7,:);
  uP = nodal_displacement(8,:);
  uQ = nodal_displacement(9,:);
  uR = nodal_displacement(10,:);

  % Construct interpolation using shape function from ANSYS Mechanical Theory Ref (and Zienkiewics Ch.4)
  interpolation_val = uI.*(2.*L1-1).*L1 + uJ.*(2.*L2-1).*L2 + uK.*(2.*L3-1).*L3 + uL.*(2.*L4-1).*L4 + ...
                      4.*uM.*L1*L2 + 4.*uN.*L2.*L3 + 4.*uO.*L1.*L3 + 4.*uP.*L1.*L4 + 4.*uQ.*L2.*L4 + ...
                      4.*uR.*L3.*L4;

end

