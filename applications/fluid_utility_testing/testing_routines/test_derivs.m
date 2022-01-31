
close all;
clear all;
clc;

addpath('../../../fluid_function_utils');

%  % Define unit tetrahedron locations
%  xI = [0,0,0];
%  xJ = [1,0,0];
%  xK = [0,1,0];
%  xL = [0,0,1];
%  xM = [0.5,0,0];
%  xN = [0.5,0.5,0];
%  xO = [0 0.5 0];
%  xP = [0 0 0.5];
%  xQ = [0.5 0 0.5];
%  xR = [0 0.5 0.5];
%  
%  nodelocs = [xI; xJ; xK; xL; xM; xN; xO; xP; xQ; xR];
%  
%  figure; hold on;
%  for m = 1:10
%    scatter3(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), 'b');
%    text(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), num2str(m));
%  end
%  xlabel('x'); ylabel('y'); zlabel('z');
%  loc1_bary = [1/3, 1/3, 1/3, 0];
%  loc2_bary = [0.6, 0.2, 0.2, 0];
%  loc3_bary = [0.2, 0.6, 0.2, 0];
%  loc4_bary = [0.2, 0.2, 0.6, 0];
%  loc1 = tet_interp(nodelocs, [1/3,1/3,1/3,0]);
%  loc2 = tet_interp(nodelocs, [0.6,0.2,0.2,0]);
%  loc3 = tet_interp(nodelocs, [0.2,0.6,0.2,0]);
%  loc4 = tet_interp(nodelocs, [0.2,0.2,0.6,0]);
%  scatter3(loc1(1), loc1(2), loc1(3), 'r');
%  scatter3(loc2(1), loc2(2), loc2(3), 'r');
%  scatter3(loc3(1), loc3(2), loc3(3), 'r');
%  scatter3(loc4(1), loc4(2), loc4(3), 'r');
%  
%  for m = 1:10
%    x = nodelocs(m,1);
%    def(m,:) = [x^2 0 0];
%  end
%  
%  nodelocs_def = nodelocs + def;
%  figure; hold on;
%  for m = 1:10
%    scatter3(nodelocs_def(m,1), nodelocs_def(m,2), nodelocs_def(m,3), 'g');
%  end
%  xlabel('x'); ylabel('y'); zlabel('z');
%  grid on;
%  
%  for m = 1:20
%    for j = 1:20
%      i = 0.5*rand(1);
%      j = 0.5*rand(1);
%      k = 1-i-j;
%      loc = tet_interp(nodelocs_def, [0,i,j,k]);
%      scatter3(loc(1), loc(2), loc(3), 'r');
%      n = determine_normal(nodelocs_def, loc, [0,i,j,k]);
%      a = nodelocs_def(2,:).*(4*i-1) + nodelocs_def(4,:).*(-4*(1-i-j)+1) + nodelocs_def(6,:).*(4*j) + nodelocs_def(9,:).*(4*(1-i-j)-4*i) + nodelocs_def(10,:).*(-4*j);
%      b = nodelocs_def(3,:).*(4*j-1) + nodelocs_def(4,:).*(-4*(1-i-j)+1) + nodelocs_def(6,:).*(4*i) + nodelocs_def(9,:).*(-4*i) + nodelocs_def(10,:).*(4*(1-i-j)-4*j);
%      nn = cross(a', b'); nn = nn ./ norm(nn);
%      disp(norm(n - nn'));
%    end
%  end
%  for m = 1:100
%    i = rand(1); j = 1-i;
%    loc = tet_interp(nodelocs_def, [0, i, j, 0]);
%    scatter3(loc(1), loc(2), loc(3), 'b');
%  end
%  xlim([0 3]); ylim([0 3]); zlim([0 3]);
%  
%  [n1] = determine_normal(nodelocs, loc1, loc1_bary);
%  [n2] = determine_normal(nodelocs, loc2, loc2_bary);
%  [n3] = determine_normal(nodelocs, loc3, loc3_bary);
%  [n4] = determine_normal(nodelocs, loc4, loc4_bary);


nodelocs = [0 0 0;
            1 0 0;
            0 1 0;
            0 0 1;
            0.5 -0.05 0;
            0 0.5 0;
            0.5 0.5 0;
            0 0 0.5;
            0.5 0 0.5;
            0 0.5 0.5];

figure; hold on; grid on;
for m = 1:10
  scatter3(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), 'b');
  text(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), num2str(m));
end
xlabel('x'); ylabel('y'); zlabel('z'); view(3);

L1list = 0.01:0.01:0.99 - 0.1;
L2list = 1 - L1list;
for m = 1:length(L1list)
  L1 = L1list(m); L2 = L2list(m);
  loc = tet_interp(nodelocs, [L1 L2 0.1 0]);
  scatter3(loc(1), loc(2), loc(3), 'r');
end

%for m = 1:10
%  x = nodelocs(m,1);
%  y = nodelocs(m,2);
%  w = 1 - (x^2 + y^2);
%  def(m,:) = [0 0 w];
%end
%nodelocs_def = nodelocs + def;
%
%loc1 = tet_interp(nodelocs, [0 1/3 1/3 1/3]);
%loc2 = tet_interp(nodelocs, [0 0.2 0.2 0.6]);
%loc3 = tet_interp(nodelocs, [0 0.2 0.6 0.2]);
%loc4 = tet_interp(nodelocs, [0 0.6 0.2 0.2]);
%
%loc1_def = tet_interp(nodelocs_def, [0 1/3 1/3 1/3]);
%loc2_def = tet_interp(nodelocs_def, [0 0.2 0.2 0.6]);
%loc3_def = tet_interp(nodelocs_def, [0 0.2 0.6 0.2]);
%loc4_def = tet_interp(nodelocs_def, [0 0.6 0.2 0.2]);
%
%figure; hold on;
%for m = 1:10
%  scatter3(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), 'b');
%  text(nodelocs(m,1), nodelocs(m,2), nodelocs(m,3), num2str(m));
%end
%xlabel('x'); ylabel('y'); zlabel('z');
%scatter3(loc1(1), loc1(2), loc1(3), 'r');
%scatter3(loc2(1), loc2(2), loc2(3), 'r');
%scatter3(loc3(1), loc3(2), loc3(3), 'r');
%scatter3(loc4(1), loc4(2), loc4(3), 'r');
%
%figure; hold on;
%for m = 1:10
%  scatter3(nodelocs_def(m,1), nodelocs_def(m,2), nodelocs_def(m,3), 'b');
%  text(nodelocs_def(m,1), nodelocs_def(m,2), nodelocs_def(m,3), num2str(m));
%end
%xlabel('x'); ylabel('y'); zlabel('z');
%scatter3(loc1_def(1), loc1_def(2), loc1_def(3), 'k');
%scatter3(loc2_def(1), loc2_def(2), loc2_def(3), 'k');
%scatter3(loc3_def(1), loc3_def(2), loc3_def(3), 'k');
%scatter3(loc4_def(1), loc4_def(2), loc4_def(3), 'k');
%
%for l = 1:100
%  i = 0.7*rand(1);
%  j = 0.3*rand(1);
%  k = 1-i-j;
%  loc = tet_interp(nodelocs_def, [0 i j k]);
%  scatter3(loc(1), loc(2), loc(3), 'k');
%  nn = determine_normal(nodelocs_def, loc, [0 i j k]);
%  r = sqrt(loc(1)^2 + loc(2)^2); t = atan2(loc(2), loc(1));
%  dwdr = -2*r;
%  na = [-dwdr*cos(t); -dwdr*sin(t); 1]; na = na ./ norm(na);
%  diff = norm(nn' - na);
%  disp(diff);
%end



function [loc] = tet_interp(nl, L)

  L1 = L(1);
  L2 = L(2);
  L3 = L(3);
  L4 = L(4);

  loc = nl(1,:).*(2*L1-1)*L1 + nl(2,:).*(2*L2-1)*L2 + nl(3,:).*(2*L3-1)*L3 + nl(4,:).*(2*L4-1)*L4 + ...
        nl(5,:).*(4*L1*L2) + nl(6,:).*(4*L2*L3) + nl(7,:).*(4*L1*L3) + nl(8,:).*(4*L1*L4) + ...
        nl(9,:).*(4*L2*L4) + nl(10,:).*(4*L3*L4);

end

