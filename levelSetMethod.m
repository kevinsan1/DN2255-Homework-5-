%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 200;
L = 2;
dx = L/(n-1);
dy = .1*dx;
% Make grid of x and y values
[boundaryX, boundaryY] = meshgrid(-L/2:L/n:L/2,-L/2:L/n:L/2);
% Define velocity field
% uext = (u,v)
% u = -cos(pi*(boundaryX+1/2)) .* sin(3*pi/8*boundaryY);
% v =  sin(pi*(boundaryX+1/2)) .* cos(3*pi/8*boundaryY);
u=2+cos(2*pi*boundaryY);
v=2+sin(2*pi*boundaryX);
tFinal = 2;
dt = 0.01;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = zeros(n+1,n+1,tSteps+1);
phi(:,:,1) = signedDistance;
phiNew = zeros(n+3,n+3,tSteps+1);
phiNew(2:(end-1),2:(end-1),1) = phi(:,:,1);
phiNew(1,:,1) = phiNew(2,:,1);
phiNew(:,1,1) = phiNew(:,2,1);
phiNew(end,:,1) = phiNew(end-1,:,1);
phiNew(:,end,1)=phiNew(:,end-1,1);
phi = phiNew;
clear phiNew;
%% Define x and y on T
theta = 0 : 2*pi/n : 2*pi;
xc = 0;
yc = -0.6;
x = r.*cos(theta) + xc;
y = r.*sin(theta) + yc;
%% Boundary Conditions
g = 2:(n+2);
ip = g + 1;
im = g - 1;
%% Plot initial phi
% figure(1)
% mesh(boundaryX,boundaryY, phi(:,:,1))
% xlabel('x')
% ylabel('y')
% i = 1;
%% Set u and v
imvel = [n+1,[1:n]];
ipvel = [[2:n+1],1];
up = u(imvel,:);
um = u(ipvel,:);
up(up<=0) = 0;
um(um>=0) = 0;
vp = v(:,imvel);vm = v(:,ipvel);
vp(vp<=0) = 0;
vm(vm>=0) = 0;
%% Main loop
for i = 1:tSteps
    wxm = phi(g,g,i) - 	phi(im,g,i);	% x backward difference
    wxp = phi(ip,g,i) - phi(g,g,i); 	% x forward difference
    wym = phi(g,g,i) - 	phi(g,im,i); 	% y backward difference
    wyp = phi(g,ip,i) - phi(g,g,i); 	% y forward difference
    phi(g,g,i+1) = phi(g,g,i)...
        -dt/dx * ( up*wxm + um*wxp )...
        -dt/dy * ( vp*wym + vm*wyp );
    phi(1,:,i+1) 		= phi(2,:,i+1);
    phi(:,1,i+1) 		= phi(:,2,i+1);
    phi(end,:,i+1) 	= phi(end-1,:,i+1);
    phi(:,end,i+1)	=phi(:,end-1,i+1);
    figure(3);clf;
    mesh(boundaryX,boundaryY, phi(g,g,i+1))
    title(sprintf('t = %0.3g',i*dt));
    xlabel('x');
    ylabel('y');
    pause(.3)
    if max(max(phi(:,:,i+1)))>1e2
        break;
    end
end
%% Plot final phi
% figure(2)
% for
% mesh(boundaryX,boundaryY, phi(:,:,i))
% xlabel('x')
% ylabel('y')
