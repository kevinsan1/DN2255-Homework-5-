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
u = setBoundariesOfNxN(u,1);
v = setBoundariesOfNxN(v,1);
tFinal = 1;
dt = 0.1*dx;
tSteps = 10;%ceil(tFinal/dt);
%% Signed Distance Function
phi = signedDistance;
phi = setBoundariesOfNxN(phi,1);
phin = zeros(n,n);

%% Set u and v
im = [1,[1:n]];
im = setBoundariesOfNxN(im,0);
ip = [[2:n],n];
ip = setBoundariesOfNxN(ip,0);
%%
up = u(im,:);
um = u(ip,:);
up(up<=0) = 0;
um(um>=0) = 0;
vp = v(:,im);vm = v(:,ip);
vp(vp<=0) = 0;
vm(vm>=0) = 0;
phiPlot(:,:,1) = phi;
%% Main loop
for tn = 1:tSteps
    for i = g
        for j = g
            wxm = phi(i,j,tn) - 	phi(im(i),j,tn);	% x backward difference
            wxp = phi(ip(i),j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - 	phi(i,im(j),tn); 	% y backward difference
            wyp = phi(i,ip(j),tn) - phi(i,j,tn); 	% y forward difference
            phin(1:n,1:n) = phi(i,j,tn)...
                -dt/dx * ( max(u(im(i),j),0)*wxm + min(u(ip(i),j),0)*wxp )...
                -dt/dy * ( max(v(i,ip(j)),0)*wym + min(v(i,ip(j)),0)*wyp );
        end
    end
    phin = setBoundariesOfNxN(phin,1);
    phiPlot(:,:,tn+1) = phin;
    
end
toc
%% Plotting
figure(2)
for iPlot = 1:tSteps
    contour(X,Y,phiPlot(g,g,1),[0,0],'r');
    axis([-1 1 -1 1])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
    pause(.001);
end

