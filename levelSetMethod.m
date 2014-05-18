%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
tic
clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 200;
dx = 2/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-1:dx:1,-1:dy:1);
g = 2:n+1; % ghost cell index
%% Define velocity field
u=-cos(pi*(X+1/2))*sin(3*pi/8*Y);
v=sin(pi*(X+1/2))*cos(3*pi/8*Y);
% u=2+cos(2*pi*Y);
% v=2+sin(2*pi*X);
tFinal = .1;
dt = .01*dx;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = signedDistance;
% phin = zeros(n+2,n+2);

%% Add Ghost Cells and Set Values
% g = 2:n+1 - ghost cell index
rbc = 'reflective velocity BC';
pbc = 'periodic BC';
nfbc = 'no flux BC';
velBC = rbc;
phiBC = pbc;
u = setBoundariesOfNxN(u,velBC);
v = setBoundariesOfNxN(v,velBC);
up(g,g) = u(g-1,g);
um(g,g) = u(g+1,g);
vp(g,g) = v(g,g-1);
vm(g,g) = v(g,g+1);
phi = setBoundariesOfNxN(phi,phiBC);
%% Set u and v
% im = [n,[1:n-1]];
% ip = [[2:n],1];
% up = u(im,:);
% um = u(ip,:);
up(up<=0) = 0;
um(um>=0) = 0;
% vp = v(:,im);vm = v(:,ip);
vp(vp<=0) = 0;
vm(vm>=0) = 0;
%% Main loop
for tn = 1:tSteps
    for i = 2:n+1
        for j = 2:n+1
            wxm = phi(i,j,tn) - phi(i-1,j,tn);	% x backward difference
            wxp = phi(i+1,j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - phi(i,j-1,tn); 	% y backward difference
            wyp = phi(i,j+1,tn) - phi(i,j,tn); 	% y forward difference
            phin(i,j) = phi(i,j,tn)...
                -dt/dx * ( up(i,j)*wxm + um(i,j)*wxp )...
                -dt/dy * ( vp(i,j)*wym + vm(i,j)*wyp );
        end
    end
    phin(1,:) = phin(2,:);
    phin(:,1) = phin(:,2);
    phin(n+2,:) = phin(n+2-1,:);
    phin(:,n+2) = phin(:,n+2-1);
    phi(:,:,tn+1) = phin(:,:);clear phin;
end
toc
%% Plotting
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
u=-cos(pi*(x+1/2))*sin(3*pi/8*y);
v=sin(pi*(x+1/2))*cos(3*pi/8*y);
for iPlot = 1:tSteps
    figure(2)
    contour(X,Y,phi(g,g,iPlot),[0,0],'r');
    hold on;
    quiver(x,y,u,v)
    hold off
    axis([-1 1 -1 1])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
    pause(.0005);
end
%%

