function [eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz] = strain(L, Jinv, q)

  % Construct B matrices
  dN1 = Jinv*L(:,1);
%  dN1 = L(:,1);
  dN1dx = dN1(1); dN1dy = dN1(2); dN1dz = dN1(3);
  B1 = [dN1dx 0 0; 0 dN1dy 0; 0 0 dN1dz; dN1dy dN1dx 0; 0 dN1dz dN1dy; dN1dz 0 dN1dx];

  dN2 = Jinv*L(:,2);
%  dN2 = L(:,2);
  dN2dx = dN2(1); dN2dy = dN2(2); dN2dz = dN2(3);
  B2 = [dN2dx 0 0; 0 dN2dy 0; 0 0 dN2dz; dN2dy dN2dx 0; 0 dN2dz dN2dy; dN2dz 0 dN2dx];

  dN3 = Jinv*L(:,3);
%  dN3 = L(:,3);
  dN3dx = dN3(1); dN3dy = dN3(2); dN3dz = dN3(3);
  B3 = [dN3dx 0 0; 0 dN3dy 0; 0 0 dN3dz; dN3dy dN3dx 0; 0 dN3dz dN3dy; dN3dz 0 dN3dx];

  dN4 = Jinv*L(:,4);
%  dN4 = L(:,4);
  dN4dx = dN4(1); dN4dy = dN4(2); dN4dz = dN4(3);
  B4 = [dN4dx 0 0; 0 dN4dy 0; 0 0 dN4dz; dN4dy dN4dx 0; 0 dN4dz dN4dy; dN4dz 0 dN4dx];

  dN5 = Jinv*L(:,5);
%  dN5 = L(:,5);
  dN5dx = dN5(1); dN5dy = dN5(2); dN5dz = dN5(3);
  B5 = [dN5dx 0 0; 0 dN5dy 0; 0 0 dN5dz; dN5dy dN5dx 0; 0 dN5dz dN5dy; dN5dz 0 dN5dx];

  dN6 = Jinv*L(:,6);
%  dN6 = L(:,6);
  dN6dx = dN6(1); dN6dy = dN6(2); dN6dz = dN6(3);
  B6 = [dN6dx 0 0; 0 dN6dy 0; 0 0 dN6dz; dN6dy dN6dx 0; 0 dN6dz dN6dy; dN6dz 0 dN6dx];

  dN7 = Jinv*L(:,7);
%  dN7 = L(:,7);
  dN7dx = dN7(1); dN7dy = dN7(2); dN7dz = dN7(3);
  B7 = [dN7dx 0 0; 0 dN7dy 0; 0 0 dN7dz; dN7dy dN7dx 0; 0 dN7dz dN7dy; dN7dz 0 dN7dx];

  dN8 = Jinv*L(:,8);
%  dN8 = L(:,8);
  dN8dx = dN8(1); dN8dy = dN8(2); dN8dz = dN8(3);
  B8 = [dN8dx 0 0; 0 dN8dy 0; 0 0 dN8dz; dN8dy dN8dx 0; 0 dN8dz dN8dy; dN8dz 0 dN8dx];

  dN9 = Jinv*L(:,9);
%  dN9 = L(:,9);
  dN9dx = dN9(1); dN9dy = dN9(2); dN9dz = dN9(3);
  B9 = [dN9dx 0 0; 0 dN9dy 0; 0 0 dN9dz; dN9dy dN9dx 0; 0 dN9dz dN9dy; dN9dz 0 dN9dx];

  dN10 = Jinv*L(:,10);
%  dN10 = L(:,10);
  dN10dx = dN10(1); dN10dy = dN10(2); dN10dz = dN10(3);
  B10 = [dN10dx 0 0; 0 dN10dy 0; 0 0 dN10dz; dN10dy dN10dx 0; 0 dN10dz dN10dy; dN10dz 0 dN10dx];

  % Assemble B matrices
  B = [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10];

  % Construct strains
  eps = B*q;

  % Unpackage strains
  eps_xx = eps(1);
  eps_yy = eps(2);
  eps_zz = eps(3);
  eps_xy = 0.5*eps(4);
  eps_yz = 0.5*eps(5);
  eps_xz = 0.5*eps(6);

end
