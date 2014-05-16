function plotCircle(xc,yc)
% plotCircle = (x,y,r) 
% (x,y) = Center Point
% r = radius
global n; global r;
theta = 0 : 2*pi/n : 2*pi;
x = r.*cos(theta) + xc;
y = r.*sin(theta) + yc;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
hold on;
plot(x,y)
hold on;
plot(xc,yc,'.')
axis([-1 1 -1 1])
% axis('square');


end

