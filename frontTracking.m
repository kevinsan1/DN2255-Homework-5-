clear all;clc;
n = 50;
grid = -1:2/(n-1):1;
[x,y] = meshgrid(grid,grid);
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
a = .3;
b = 1;
xplus1 = x(:,ip);
xmin1 = x(:,im);
yplus1 = y(ip,:);
ymin1 = y(im,:);
xden = xplus1-x;
yden = yplus1-y;
uabc = 1./(xden.*yden');
for i = 1:n
    for j = 1:n
        uab(i,j) = ...
            u(i,j)*(xv(ip(i))-a)*(yv(jp(j))-b) + ...
            u(ip(i),j)* (a-xv(i)) * (yv(jp(j)) - b) + ...
            u(i,jp(j))*(xv(ip(i))-a)*(b-yv(j)) + ...
            u(ip(i),jp(j))*(a-xv(i))*(b-yv(j));
    end
end
