clearvars;clc;close all;
load('t1p6.mat');
num = 20;
[C,h] = contour(X,Y,phi(:,:,num),[0,0],'b');
%%
% hold on;
x = C(1,2:end);
y = C(2,2:end);
n = length(x);
L=2;
%% Make grid of x and y values
[boundaryX, boundaryY] = meshgrid(-L/2:L/(n-1):L/2,-L/2:L/(n-1):L/2);
%% Preallocation
d = zeros(1,n);
minDist = zeros(n,n);
%% Get minimum distance from (boundaryX,boundaryY) to any (x,y)
% for i = 1:(n)
%     for k = 1:(n)
%         for j = 1:(n)
%             d(j) = (boundaryX(i,k) - x(j)).^2 + ...
%                 (boundaryY(i,k) - y(j)).^2;
%         end
%         [minDist(i,k),index] = min(sqrt(d));
%     end
% end
% %% Test if inside circle
% for i = 1:n
%     for j = 1:n
%         if testnumber(i,j) <= r^2
%             corrMinDist(i,j) = minDist(i,j);
%             insideCircleTest(i,j) = 1;
%         elseif testnumber(i,j) > r^2
%             corrMinDist(i,j) = -1*minDist(i,j);
%         end
%     end
% end
xcon = -L/2:L/(n-1):L/2;
ycon = xcon;
% contour(xcon,ycon,phiCorrect,[0,0]);
% hold off
% legend('\phi computed','phiCorrect');
%%
% figure(2);clf
% surf(X,Y,phi(:,:,num))
axis([-1 1 0.2 1 -0.02 0.02]);
xlabel('x');
ylabel('y');
phin = phi(:,:,num);
phin(phi(:,:,num)<=0)=0;
isequal(phin,zeros(length(phin)))
%%
contour(X,Y,phin,[0,0],'b')
hold on;
contour(X,Y,phi(:,:,num),[0,0],'r');
% contour(xcon,ycon,phiCorrect,[0,0],'k');
legend('\phi computed','phiCorrect');
hold off;
%%
figure(4);clf
mesh(phin)
hold on;
surf(phi(:,:,num));
hold off;