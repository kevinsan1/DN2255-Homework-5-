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
c = .01;
dt = 0.001*dx/c;
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
for t = 1:ceil(tfinal/dt)
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
        xTn(istep) = xT(istep) - c*dt/dx*uab;
        yTn(istep) = yT(istep) - c*dt/dx*vab;
    end
    xT = xTn;
    yT = yTn;
    xplot(:,t+1) = xT;
    yplot(:,t+1) = yT;
end
%% matrix way
for tm = 1:10
    for i = 1:n
        a(i)=find(x(1,:)<=xT(i),1,'last');
        b(i)=find(y(:,1)<=yT(i),1,'last');
    end
    figure(1);clf
    plot(x(1,a),y(b,1),'.','color','r')
    hold on;
    plot(xT,yT,'.')
    axis([-.4,.4,0.2,1])
    axis('square')
    hold off;
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
    xTn(istep) = xT(istep) - c*dt/dx*uab;
    yTn(istep) = yT(istep) - c*dt/dx*vab;
end
%%
% xplus1 = x(:,ip);
% xmin1 = x(:,im);
% yplus1 = y(ip,:);
% ymin1 = y(im,:);
% xden = xplus1-x;
% yden = yplus1-y;
% uabc = 1./(xden.*yden');
%% Plot
q = 1:n/20:n;
xq = x(q,q);
yq = y(q,q);
uq =  u(q,q);
vq =  v(q,q);
figure(1);clf;
for ix = 1:t-1
    plot(xplot(:,ix),yplot(:,ix))
    hold on;
    quiver(xq,yq,uq,vq)
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    hold off;
    pause(0.001)
end


