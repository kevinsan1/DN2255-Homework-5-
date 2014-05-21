% ix = find(xv<=a,1,'last');
% xv(ix)
% xplus1 = x(:,ip);
% xmin1 = x(:,im);
% yplus1 = y(ip,:);
% ymin1 = y(im,:);
% xden = xplus1-x;
% yden = yplus1-y;
% uabc = 1./(xden.*yden');
theta = 0 : 2*pi/(n-1) : 2*pi;
xT = r.*cos(theta) + xc;
yT = r.*sin(theta) + yc;
figure(1);clf;
hold on;
axis([-1,1,-1,1])
axis('square')
for i = 1:n
    plot(xT(i),yT(i),'.')
    pause(0.01)
end