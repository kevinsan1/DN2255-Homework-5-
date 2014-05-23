clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
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
% q = 1:n/20:n;
% xq = x(q,q);
% yq = y(q,q);
% uq =  u(q,q);
% vq =  v(q,q);
% figure(1);clf;
% for ix = 1:t-1
%     plot(oxplot(:,ix),oyplot(:,ix),'.',xplot(:,ix),yplot(:,ix),'.')
%     hold on;
%     quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix));
%     quiver(xq,yq,uq,vq)
%     axis([-L/2 L/2 -L/2 L/2])
%     axis('square')
%     hold off;
%     pause(0.001)
% end
%% Plot of discrete times figure 2
quivLinewidth = 1.1;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
ix = 1;
p1=plot(oxplot(:,ix),oyplot(:,ix),'.')
hold on;
q1 = quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix)); % quivPlots
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
hTitle1=title(sprintf('$t$ = %g',dt*ix-dt));
hLegend1 = legend('($x_{i}, y_{j})$','$(u_{ab},v_{ab})$','$(u,v)$','location','best');
axis('square')
set(p1,'MarkerSize',10)
set(q1,'LineWidth',quivLinewidth);
hold off;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
ix = 25;
p2=plot(oxplot(:,ix),oyplot(:,ix),'.')
hold on;
q2 = quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix)); % quivPlots
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
hTitle2 = title(sprintf('$t$ = %g',dt*ix-dt));
hLegend2 = legend('($x_{i}, y_{j})$','$(u_{ab},v_{ab})$','$(u,v)$','location','best');
axis('square')
set(p2,'MarkerSize',10)
set(q2,'LineWidth',quivLinewidth);
hold off;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
ix = 50;
p3=plot(oxplot(:,ix),oyplot(:,ix),'.')
hold on;
q3 = quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix)); % quivPlots
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
hTitle3 = title(sprintf('$t$ = %g',dt*ix-dt));
hLegend3 = legend('($x_{i}, y_{j})$','$(u_{ab},v_{ab})$','$(u,v)$','location','best');
axis('square')
set(p3,'MarkerSize',10)
set(q3,'LineWidth',quivLinewidth);
hold off;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
ix = 100;
p4=plot(oxplot(:,ix),oyplot(:,ix),'.')
hold on;
q4 = quiver(oxplot(:,ix),oyplot(:,ix),quplot(:,ix),qvplot(:,ix)); % quivPlots
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
hTitle4 = title(sprintf('$t$ = %g',dt*ix-dt));
hLegend4 = legend('($x_{i}, y_{j})$','$(u_{ab},v_{ab})$','$(u,v)$','location','best');
set(p4,'MarkerSize',10)
set(q4,'LineWidth',quivLinewidth);
set([hTitle1, hTitle2, hTitle3,hTitle4,...
    hLegend1, hLegend2, hLegend3,hLegend4],...
    	'FontSize'   , 14          , ...
    	'FontWeight' , 'bold'      ,...
        'FontName'   , 'Helvetica',...
        'Interpreter','Latex');
axis('square')
hold off;
printYesNo = 1;
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
if printYesNo == 1
    figure(1);
    print('-depsc2', [saveFigurePath ...
        sprintf('fronttrackingfig1')]);
    figure(2)
    print('-depsc2', [saveFigurePath ...
        sprintf('fronttrackingfig2')]);
    figure(3)
    print('-depsc2', [saveFigurePath ...
        sprintf('fronttrackingfig3')]);
    figure(4)
    print('-depsc2', [saveFigurePath ...
        sprintf('fronttrackingfig4')]);
    makeTable({'dt';'dx';'n'},[dt, dx, n]','myTablefronttracking.tex')
end


