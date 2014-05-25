function corrMinDist = testIfInsideCircle(x,y,minDist)
%  	Description
%	corrMinDist = testIfInsideCircle(x)
%   x = grid points
%	y = grid points
%   circleX = x center of circle
%   circleY = y center of circle
%	minDist = min distance from x to xc and y to yc
global n r insideCircleTest yc xc;

testnumber = (x-xc).^2 + (y-yc).^2;
% corrMinDist = minDist;
% pointsInside = cell(n+1);
for i = 1:n
    for j = 1:n
        if testnumber(i,j) <= r^2
            corrMinDist(i,j) = minDist(i,j);
            insideCircleTest(i,j) = 1;
        elseif testnumber(i,j) > r^2
            corrMinDist(i,j) = -1*minDist(i,j);
        end
    end
end

