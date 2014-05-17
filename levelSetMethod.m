%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
tic
clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 200;
dx = 2/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-1:dx:1,-1:dy:1);
% Define velocity field
% uext = (u,v)
% u = -cos(pi*(boundaryX+1/2)) .* sin(3*pi/8*boundaryY);
% v =  sin(pi*(boundaryX+1/2)) .* cos(3*pi/8*boundaryY);

g = 2:n+1;
u=2+cos(2*pi*Y);
v=2+sin(2*pi*X);
u = setBoundariesOfNxN(u);
v = setBoundariesOfNxN(v);
tFinal = 1;
dt = 0.1*dx;
tSteps = 10;%ceil(tFinal/dt);
%% Signed Distance Function
phi = signedDistance;
phi = setBoundariesOfNxN(phi);
phin = zeros(n,n);

%% Set u and v
im = [1,[1:n],n];
ip = [2,[2:n],n];
up = u(im,:);
um = u(ip,:);
up(up<=0) = 0;
um(um>=0) = 0;
vp = v(:,im);vm = v(:,ip);
vp(vp<=0) = 0;
vm(vm>=0) = 0;
%% Main loop
for tn = 1:tSteps
    wxm = phi(g,g,tn) - 	phi(g-1,g,tn);	% x backward difference
    wxp = phi(g+1,g,tn) - phi(g,g,tn); 	% x forward difference
    wym = phi(g,g,tn) - 	phi(g,g-1,tn); 	% y backward difference
    wyp = phi(g,g+1,tn) - phi(g,g,tn); 	% y forward difference
    phin = phi(g,g,tn)...
        -dt/dx * ( max(u(g,g),0)*wxm + min(u(g,g),0)*wxp )...
        -dt/dy * ( max(v(g,g),0)*wym + min(v(g,g),0)*wyp );
    phin = setBoundariesOfNxN(phin);
    phi(:,:,tn+1) = phin;
end
toc
%% Plotting
figure(2)
for iPlot = 1:tSteps
    contour(X,Y,phi(g,g,iPlot),[0,0],'r');
    axis([-1 1 -1 1])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
    pause(.001);
end

