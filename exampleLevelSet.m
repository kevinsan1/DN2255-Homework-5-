%% Two dimensional level method for a prescribed
% velocity field based on upwinding
%
clear all
close all
N = 100;         % number of grid points in one direction
numc = 1;        % number of circles for the initial condition
R1 = .3;         % initial radius of circles
h = 2/(N-1);     % grid spacing
dt = .2*h;       % time step
tfin = 1.6;       % total simulation time
nit = tfin/dt;    % number of time steps
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
u=-cos(pi*(X+.5)).*sin(3*pi/8*Y);
v=sin(pi*(X+.5)).*cos(3*pi/8*Y);
qp = 1:5:N;
%      arrays for the periodic boundary conditions
for i=1:N
    ip(i)=i+1;
    im(i)=i-1;
end
im(1)=N;
ip(N)=1;
%
% begin simulation loop
for iter=1:nit
    for i=1:N
        for j=1:N
            dmx=(phi(i,j)-phi(im(i),j))/h; % x backward difference
            dpx=(phi(ip(i),j)-phi(i,j))/h; % x forward difference
            dmy=(phi(i,j)-phi(i,im(j)))/h; % y backward difference
            dpy=(phi(i,ip(j))-phi(i,j))/h; % y forward difference
            convx=max(u(i,j),0)*dmx+min(u(i,j),0)*dpx;
            convy=max(v(i,j),0)*dmy+min(v(i,j),0)*dpy;
            phin(i,j)=phi(i,j)-(convx+convy)*dt; % advance by dt
        end
    end
    phi=phin; % update
    %
    % Plotting
    contour(X,Y,phi,[0,0],'r');
    hold on;
    quiver(X(qp,qp),Y(qp,qp),u(qp,qp),v(qp,qp))
    title(sprintf('t = %g', iter*dt))
    axis([-1 1 -1 1])
    axis('square')
    hold off;
    pause(.001)
end