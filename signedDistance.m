clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/'];
addpath(genpath(myPath));
%% Define constants
n = 200;
%% Define x and y on T
rMax = 0.3;
theta = 0 : 2*pi/n : 2*pi;
r = rMax;
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

%% Plot minDist
% mesh(boundaryX, boundaryY,minDist)
% hold on;
% mesh(x,y,zeros(n+1,n+1))
% %% Plot initial circle
% plotCircle(xc,yc,r,n);
