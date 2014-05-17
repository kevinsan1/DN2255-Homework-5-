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
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:L/n:L/2,-L/2:L/n:L/2);
% Define velocity field
% uext = (u,v)
% u = -cos(pi*(boundaryX+1/2)) .* sin(3*pi/8*boundaryY);
% v =  sin(pi*(boundaryX+1/2)) .* cos(3*pi/8*boundaryY);
u=2+cos(2*pi*Y);
v=2+sin(2*pi*X);
tFinal = 2;
dt = 0.1*dx;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = zeros(n+1,n+1,tSteps+1);
phi(:,:,1) = signedDistance;
phin = zeros(n+3,n+3,tSteps+1);
phin(2:(end-1),2:(end-1),1) = phi(:,:,1);
phin(1,:,1) = phin(2,:,1);
phin(:,1,1) = phin(:,2,1);
phin(end,:,1) = phin(end-1,:,1);
phin(:,end,1)=phin(:,end-1,1);
phi = phin;
phin = zeros(n+3,n+3);
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
%         Plotting
contour(X,Y,phi(g,g,1),[0,0],'r');
axis([-1 1 -1 1])
axis('square')
pause(.001)
tic
%% Main loop
for tn = 1:tSteps
    for i = 2:n+2
        for j = 2:n+2
            wxm = phi(i,j,tn) - 	phi(imvel(i-1),j,tn);	% x backward difference
            wxp = phi(ipvel(i-1),j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - 	phi(i,imvel(j-1),tn); 	% y backward difference
            wyp = phi(i,ipvel(j-1),tn) - phi(i,j,tn); 	% y forward difference
            phin(i,j) = phi(i,j,tn)...
                -dt/dx * ( max(u(i-1,j-1),0)*wxm + min(u(i-1,j-1),0)*wxp )...
                -dt/dy * ( max(v(i-1,j-1),0)*wym + min(v(i-1,j-1),0)*wyp );
            %             figure(3);clf;
            %             mesh(boundaryX,boundaryY, phi(i,j,tn+1))
            %             title(sprintf('t = %0.3g',tn*dt));
            %             xlabel('x');
            %             ylabel('y');
            %             pause(.3)
            if max(max(phi(:,:,tn+1)))>1e2
                break;
            end
        end
    end
    phin(1,:) 		= phin(2,:);
    phin(:,1) 		= phin(:,2);
    phin(end,:) 	= phin(end-1,:);
    phin(:,end)     = phin(:,end-1);
    phi(:,:,tn+1) = phin;
    figure(2)
    contour(X,Y,phi(g,g,tn),[0,0],'r');
    axis([-1 1 -1 1])
    axis('square')
    pause(.001);
end
toc
%% Plot final phi
%         Plotting
for iPlot = 1:131
    figure(2)
    contour(X,Y,phi(g,g,iPlot),[0,0],'r');
    axis([-1 1 -1 1])
    axis('square')
    pause(.001)
end
