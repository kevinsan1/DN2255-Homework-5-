clear all;clc;myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/'];
addpath(genpath(myPath));
%% Define x and y on T
thetaMax = 2*pi;
rMax = 1;
n = 100;
theta = 0 : thetaMax/n : thetaMax;
r = rMax;
x = r.*cos(theta);
y = r.*sin(theta);
%% Plot figure
% figure(1);clf
% plot(x,y)
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



