function val = localeval3d(x,y,z,nodesval)
% Evaluate values inside the rectangular element. 0 <= x,y,z <= 1 is the
% local coordinate. nodesval contains the values at the four vertex in
% counter-clockwise and down-first order:
% ======== My beautiful drawing of a cube ======
%    6-------7
%   /|      /|    |
%  / |     / |    z
% 5-------8  |    |
% |  2----|--3
% | /     | /   /
% |/      |/   y
% 1-------4   /
%  
% ----x---
% ==============================================
% xxi is on the x-axis and eeta y-axis

%% Basic error checking
if x < 0
    error('x < 0');
end
if y < 0
    error('y < 0');
end
if z < 0
    error('z < 0');
end
if x > 1
    error('x> 1');
end
if y > 1
    error('y > 1');
end
if z > 1
    error('z > 1');
end

A = ...
   [ 1     0     0     0     0     0     0     0; ...
    -1     1     0     0     0     0     0     0; ...
    -1     0     1     0     0     0     0     0; ...
     1    -1    -1     1     0     0     0     0; ...
    -1     0     0     0     1     0     0     0; ...
     1    -1     0     0    -1     1     0     0; ...
     1     0    -1     0    -1     0     1     0; ...
    -1     1     1    -1     1    -1    -1     1];
basis = [x^0*y^0*z^0 x^1*y^0*z^0 x^0*y^1*z^0 x^1*y^1*z^0 x^0*y^0*z^1 x^1*y^0*z^1 x^0*y^1*z^1 x^1*y^1*z^1];

val = basis*A*nodesval;
end