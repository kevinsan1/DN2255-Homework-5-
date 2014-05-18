[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
% u = 2+cos(2*pi*x);
% v = 2+sin(2*pi*y);
u=-cos(pi*(x+1/2))*sin(3*pi/8*y);
v=sin(pi*(x+1/2))*cos(3*pi/8*y);
figure
quiver(x,y,u,v)