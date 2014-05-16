function corrMinDist = signedDistance
%  	Description
%	corrMinDist = signedDistance
global n; global r; global insideCircleTest;
%% Define global constants
insideCircleTest = zeros(n+1);
r = 0.3;
%% Define x and y on T
theta = 0 : 2*pi/n : 2*pi;
xc = 0;
yc = -0.6;
x = r.*cos(theta) + xc;
y = r.*sin(theta) + yc;
%% Make grid of x and y values
[boundaryX, boundaryY] = meshgrid(-1:2/n:1,-1:2/n:1);
%% Preallocation
d = zeros(1,n+1);
minDist = zeros(n+1,n+1);
%% Get minimum distance from (boundaryX,boundaryY) to any (x,y)
for i = 1:(n+1)
    for k = 1:(n+1)
        for j = 1:(n+1)
            d(j) = (boundaryX(i,k) - x(j)).^2 + ...
            (boundaryY(i,k) - y(j)).^2;
        end
        [minDist(i,k),index] = min(d);
    end
end
clear i j k;
%% Set minDist to negative if inside circle T - (x,y)
corrMinDist = testIfInsideCircle(boundaryX, ...
 boundaryY, minDist, xc, yc);
%% Plot minDist
% mesh(boundaryX, boundaryY,corrMinDist)
% hold on;
% mesh(boundaryX, boundaryY,minDist)
% xlabel('x')
% ylabel('y')
% diff = corrMinDist - minDist;



end % function