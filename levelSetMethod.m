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
u(g,g)=2+cos(2*pi*Y);
v(g,g)=2+sin(2*pi*X);
u(1,:)      = u(2,:);
u(end,:)    = u(end-1,:);
v(1,:)      = v(2,:);
v(end,:)    = v(end-1,:);
tFinal = 1;
dt = 0.1*dx;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = zeros(n+2,n+2);
phi(g,g,1) = signedDistance;
phin = zeros(n+2,n+2);
%% Define x and y on T
theta = 0 : 2*pi/n : 2*pi;
xc = 0;
yc = -0.6;
x = r.*cos(theta) + xc;
y = r.*sin(theta) + yc;
%% Plot initial phi
% figure(1)
% mesh(boundaryX,boundaryY, phi(:,:,1))
% xlabel('x')
% ylabel('y')
% i = 1;
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
    for i = g
        for j = g
            wxm = phi(i,j,tn) - 	phi(im(i),j,tn);	% x backward difference
            wxp = phi(ip(i),j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - 	phi(i,im(j),tn); 	% y backward difference
            wyp = phi(i,ip(j),tn) - phi(i,j,tn); 	% y forward difference
            phin(i,j) = phi(i,j,tn)...
                -dt/dx * ( max(u(im(i),j),0)*wxm + min(u(ip(i),j),0)*wxp )...
                -dt/dy * ( max(v(i,ip(j)),0)*wym + min(v(i,ip(j)),0)*wyp );
        end
    end
    phin(1,:) 		= phin(2,:);
    phin(:,1) 		= phin(:,2);
    phin(end,:) 	= phin(end-1,:);
    phin(:,end)     = phin(:,end-1);
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

