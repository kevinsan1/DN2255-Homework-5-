function corrMinDist = signedDistance(x,y)
%  	Description
%	corrMinDist = signedDistance
global n xc yc r insideCircleTest;
%% Define global constants
% insideCircleTest = zeros(n+1);

%% Define x and y on T
if nargin < 2
    theta = 0 : 2*pi/(n-1) : 2*pi;
    x = r.*cos(theta) + xc;
    y = r.*sin(theta) + yc;
end
L = 2;
%% Make grid of x and y values
[boundaryX, boundaryY] = meshgrid(-L/2:L/(n-1):L/2,-L/2:L/(n-1):L/2);
%% Preallocation
d = zeros(1,n);
minDist = zeros(n,n);
%% Get minimum distance from (boundaryX,boundaryY) to any (x,y)
for i = 1:(n)
    for k = 1:(n)
        for j = 1:(n)
            d(j) = (boundaryX(i,k) - x(j)).^2 + ...
                (boundaryY(i,k) - y(j)).^2;
        end
        [minDist(i,k),index] = min(sqrt(d));
    end
end
clear i j k;
%% Set minDist to negative if inside circle T - (x,y)
corrMinDist = testIfInsideCircle(boundaryX, ...
    boundaryY, minDist);
%% Plot minDist
% mesh(boundaryX, boundaryY,corrMinDist)
% hold on;
% mesh(boundaryX, boundaryY,minDist)
% xlabel('x')
% ylabel('y')
% diff = corrMinDist - minDist;



end % function