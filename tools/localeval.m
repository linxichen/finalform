function val = localeval3d(x,y,z,nodesval)
% Evaluate values inside the rectangular element. 0 <= x,y,z <= 1 is the
% local coordinate. nodesval contains the values at the four vertex in
% counter-clockwise and down-first order:
% ======== My beautiful drawing of a cube ======
%    6-------7
%   /|      /|
%  / |     / |
% 5-------8  |
% |  2----|--3
% | /     | /
% |/      |/
% 1-------4
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

val1 = nodesval(1);
val2 = nodesval(2); 
val3 = nodesval(3); 
val4 = nodesval(4); 
val = (1-xxi)*(1-eeta)*val1 + xxi*(1-eeta)*val2 + xxi*eeta*val3 + (1-xxi)*eeta*val4;
end