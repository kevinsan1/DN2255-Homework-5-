clear all;clc;
global n xc yc r
%% Definitions
n = 100;
xc = 0;
yc = 0.6;
r = 0.3;
L = 2;
dx = L/(n-1);
grid = -L/2:dx:L/2;
[x,y] = meshgrid(grid,grid);
c = 1;
dt = 0.01*dx/c;
tfinal = 1.6;
%% Velocity Field
u=-cos(pi*(x+0.5)).*sin(3*pi/8*y);
v=sin(pi*(x+0.5)).*cos(3*pi/8*y);
%% Boundary Conditions
ip = [2:n,1];
jp = ip;
im = [n,1:n-1];
jm = im;
xv = grid;
yv = grid;
%% Velocity at point (a,b)
a = 10;
b = 12;
ix = find(xv<=a,1,'last');
jx = find(yv<=b,1,'last');
uab = ...
    1/((xv(ip(ix))-xv(ix))*(yv(jp(jx))-yv(jx)))*...
    (u(ix,jx)*(xv(ip(ix))-a)*(yv(jp(jx))-b) + ...
    u(ip(ix),jx)* (a-xv(ix)) * (yv(jp(jx)) - b) + ...
    u(ix,jp(jx))*(xv(ip(ix))-a)*(b-yv(jx)) + ...
    u(ip(ix),jp(jx))*(a-xv(ix))*(b-yv(jx)));
vab = ...
    1/((xv(ip(ix))-xv(ix))*(yv(jp(jx))-yv(jx)))*...
    (v(ix,jx)*(xv(ip(ix))-a)*(yv(jp(jx))-b) + ...
    v(ip(ix),jx)* (a-xv(ix)) * (yv(jp(jx)) - b) + ...
    v(ix,jp(jx))*(xv(ip(ix))-a)*(b-yv(jx)) + ...
    v(ip(ix),jp(jx))*(a-xv(ix))*(b-yv(jx)));
%% Define x and y on T
theta = 0 : 2*pi/(n-1) : 2*pi;
xT = r.*cos(theta) + xc;
yT = r.*sin(theta) + yc;
xplot(:,1) = xT;
yplot(:,1) = yT;
%% Advection
q = 1:n/20:n;
xq = x(q,q);
yq = y(q,q);
uq =  u(q,q);
vq =  v(q,q);
for t = 1:114
    for istep = 1:n
        ix = find(xv<=xT(istep),1,'last');
        jx = find(yv<=yT(istep),1,'last');
        uab = ...
            1/((xv(ip(ix))-xv(ix))*(yv(jp(jx))-yv(jx)))*...
            (u(ix,jx)*(xv(ip(ix))-xT(istep))*(yv(jp(jx))-yT(istep)) + ...
            u(ip(ix),jx)* (xT(istep)-xv(ix)) * (yv(jp(jx)) - yT(istep)) + ...
            u(ix,jp(jx))*(xv(ip(ix))-xT(istep))*(yT(istep)-yv(jx)) + ...
            u(ip(ix),jp(jx))*(xT(istep)-xv(ix))*(yT(istep)-yv(jx)));
        vab = ...
            1/((xv(ip(ix))-xv(ix))*(yv(jp(jx))-yv(jx)))*...
            (v(ix,jx)*(xv(ip(ix))-xT(istep))*(yv(jp(jx))-yT(istep)) + ...
            v(ip(ix),jx)* (xT(istep)-xv(ix)) * (yv(jp(jx)) - yT(istep)) + ...
            v(ix,jp(jx))*(xv(ip(ix))-xT(istep))*(yT(istep)-yv(jx)) + ...
            v(ip(ix),jp(jx))*(xT(istep)-xv(ix))*(yT(istep)-yv(jx)));
        xTn(istep) = xT(istep) + c*dt/dx*uab;
        yTn(istep) = yT(istep) + c*dt/dx*vab;
        oxplot(istep,t) = xv(ix);
        oyplot(istep,t) = yv(jx);
        quplot(istep,t) = uab;
        qvplot(istep,t) = vab;
    end
    xT = xTn;
    yT = yTn;
    xplot(:,t+1) = xT;
    yplot(:,t+1) = yT;
    %     figure(1);
    %     plot(oxplot,oyplot,'.',xplot(:,t),yplot(:,t),'.');
    %     hold on;
    %     quiver(oxplot,oyplot,quplot,qvplot);
    %     quiver(xq,yq,uq,vq);
    %     axis([-1,1,0,1])
    %     axis('square')
    %     hold off;
    %     pause(0.001)
end
%% Plot of all time figure 1
q = 1:n/20:n;
xq = x(q,q);
yq = y(q,q);
uq =  u(q,q);
vq =  v(q,q);
figure(1);clf;
for ix = 1:t-1
    plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'.')
    hold on;
    quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
    quiver(xq,yq,uq,vq)
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    hold off;
    pause(0.001)
end
%% Plot of discrete times figure 2
figure(2)
mn = 4;
nn = 1;
subplot(mn,nn,1) % subplot 1
ix = 1;
plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'-')
hold on;
quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
title(sprintf('t = %g',dt*ix-dt))
hLegend = legend('a','(x_{i}, y_{j})','c','d','location','best')
axis('square')
hold off;
pause(0.001)
subplot(mn,nn,2) % subplot 2
ix = 25;
plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'.')
hold on;
quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
title(sprintf('t = %g',dt*ix-dt))
axis('square')
hold off;
pause(0.001)
subplot(mn,nn,3) % subplot 3
ix = 50;
plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'.')
hold on;
quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
title(sprintf('t = %g',dt*ix-dt))
axis('square')
hold off;
pause(0.001)
subplot(mn,nn,4) % subplot 4
ix = 100;
plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'.')
hold on;
quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
title(sprintf('t = %g',dt*ix-dt))
axis('square')
hold off;
pause(0.001)


