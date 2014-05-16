%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 100;
L = 2;
dx = L/n;
dy = dx;
% Make grid of x and y values
[boundaryX, boundaryY] = meshgrid(-L/2:L/n:L/2,-L/2:L/n:L/2);
% Define velocity field
% uext = (u,v)
u = -cos(pi*(boundaryX+1/2)) .* sin(3*pi/8*boundaryY);
v =  sin(pi*(boundaryX+1/2)) .* cos(3*pi/8*boundaryY);
a = sqrt(u.^2 + v.^2);
tFinal = 1.6;
aAv = sum(sum(a))/length(a)^2;
dt = dx/aAv;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = zeros(n+1,n+1,tSteps+1);
phi(:,:,1) = signedDistance;
phi(:,:,2) = phi(:,:,1);
%% Define x and y on T
theta = 0 : 2*pi/n : 2*pi;
xc = 0;
yc = -0.6;
x = r.*cos(theta) + xc;
y = r.*sin(theta) + yc;
%% Boundary Conditions
% 
%% Plot initial phi
figure(1)
mesh(boundaryX,boundaryY, phi(:,:,1))
xlabel('x')
ylabel('y')
% i = 1;
%% Set u and v
up = u;
um = u;
up(up<=0) = 0;
um(um>=0) = 0;
vp = v;
vm = v;
vp(vp<=0) = 0;
vm(vm>=0) = 0;
%% Main loop
for i = 2:tSteps
    W(i-1) = (phi(:,:,i) - phi(:,:,i-1));
    W(i+1) = (phi(:,:,i) - phi(:,:,i+1));
    phi(:,:,i+1) = phi(:,:,i)...
        -dt/dx * ( up(im)*W(i-1) + um(ip)*W(i+1))...
        -dt/dy * ()
    figure(3);clf;
    mesh(boundaryX,boundaryY, phi(:,:,i))
    title(sprintf('t = %0.3g',i*dt));
    pause(0.3)
end
%% Plot final phi
% figure(2)
% for
% mesh(boundaryX,boundaryY, phi(:,:,i))
% xlabel('x')
% ylabel('y')
