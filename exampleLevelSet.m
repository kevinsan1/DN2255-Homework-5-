% Level Set Functions for Certain Shapes
%
% Two dimensional level method for a prescribed
% velocity field based on upwinding
%
clear all
close all
N=200;         % number of grid points in one direction
numc=1;        % number of circles for the initial condition
R1=.3;         % initial radius of circles
h=2/(N-1);     % grid spacing
dt=.1*h;       % time step
tfin=2;       % total simulation time
nit=tfin/dt;   % number of time steps
x=-1:h:1;
y=x;
[X,Y]=meshgrid(x);
%
%
%   Initialize the level set function
%   Union of numc randomly located circles
a1=2*(1-R1-h)*rand-(1-R1-h);
b1=2*(1-R1-h)*rand-(1-R1-h);
phi=((X-a1).*(X-a1)+(Y-b1).*(Y-b1)).^.5-R1;
for k=1:numc-1
    a1=2*(1-R1-h)*rand-(1-R1-h);
    b1=2*(1-R1-h)*rand-(1-R1-h);
    phi1=((X-a1).*(X-a1)+(Y-b1).*(Y-b1)).^.5-R1;
    phi=min(phi,phi1);
end
figure(1)
contour(X,Y,phi,[0,0],'r');
axis([-1 1 -1 1])
axis('square')
%
%   Initialize the velocity field
    u=2-cos(2*pi*Y);
    v=2+sin(2*pi*X);
%
%      arrays for the periodic boundary conditions
for i=1:N
    ip(i)=i+1;
    im(i)=i-1;
end
im(1)=N;
ip(N)=1;
%
%      begin simulation loop
for iter=1:nit
    for i=1:N
        for j=1:N
            dmx=(phi(i,j)-phi(im(i),j))/h;                     % x backward difference
            dpx=(phi(ip(i),j)-phi(i,j))/h;                     % x forward difference
            dmy=(phi(i,j)-phi(i,im(j)))/h;                     % y backward difference
            dpy=(phi(i,ip(j))-phi(i,j))/h;                     % y forward difference
            convx=max(u(i,j),0)*dmx+min(u(i,j),0)*dpx;
            convy=max(v(i,j),0)*dmy+min(v(i,j),0)*dpy;
            phin(i,j)=phi(i,j)-(convx+convy)*dt;         % advance by dt
        end
    end
    phi=phin;                                             % update
    %
    %         Plotting
    contour(X,Y,phi,[0,0],'r');
    axis([-1 1 -1 1])
    axis('square')
    pause(.001)
end