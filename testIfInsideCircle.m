function corrMinDist = testIfInsideCircle(x,y,minDist,circleX, circleY)
%  	Description
%	corrMinDist = testIfInsideCircle(x)
%   x = 
%	y =
%   circleX
%   circleY
%	minDist =
global n; global r; global insideCircleTest;
testnumber = (x-circleX).^2 + (y-circleY).^2;
corrMinDist = minDist;
% pointsInside = cell(n+1);
for i = 1:n+1
    for j = 1:n+1
        if testnumber(i,j) < r^2
            corrMinDist(i,j) = minDist(i,j);
            insideCircleTest(i,j) = 1;
        elseif testnumber(i,j) > r^2
            corrMinDist(i,j) = -1*minDist(i,j);
        elseif testnumber(i,j) == r^2
            corrMinDist(i,j) = minDist(i,j);
        end
    end
end
end % function
