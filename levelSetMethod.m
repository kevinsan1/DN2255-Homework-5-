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
tic
%% Main loop
for tn = 1:tSteps
    for i = 2:n+2
        for j = 2:n+2
            wxm = phi(i,j,tn) - 	phi(i-1,j,tn);	% x backward difference
            wxp = phi(i+1,j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - 	phi(i,j-1,tn); 	% y backward difference
            wyp = phi(i,j+1,tn) - phi(i,j,tn); 	% y forward difference
            phi(i,j,tn+1) = phi(i,j,tn)...
                -dt/dx * ( up(i-1,j-1)*wxm + um(i-1,j-1)*wxp )...
                -dt/dy * ( vp(i-1,j-1)*wym + vm(i-1,j-1)*wyp );
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
    phi(1,:,tn+1) 		= phi(2,:,tn+1);
    phi(:,1,tn+1) 		= phi(:,2,tn+1);
    phi(end,:,tn+1) 	= phi(end-1,:,tn+1);
    phi(:,end,tn+1)	=phi(:,end-1,tn+1);
end
toc
%% Plot final phi
figure(2)
for i = 1:tSteps
    pause(0.5)
    mesh(boundaryX,boundaryY, phi(g,g,i))
    xlabel('x')
    ylabel('y')
end
