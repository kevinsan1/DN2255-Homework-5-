clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
global n xc yc r
%% Definitions
n = 100;
xc = 0;
yc = -0.6;
r = 0.3;
L = 2;
dx = L/(n-1);
vec = -L/2:dx:L/2;
[y,x] = meshgrid(vec,vec);
c = 1;
dt = 0.1*dx;
tfinal = 1.6;
totStep = ceil(tfinal/dt);
%% Velocity Field
u=-cos(pi*(x+0.5)).*sin(3*pi/8*y);
v=sin(pi*(x+0.5)).*cos(3*pi/8*y);
xv = vec;
yv = vec;
%% Define x and y on T
theta = 0 : 2*pi/(n-1) : 2*pi;
xT = r.*cos(theta) + xc;
yT = r.*sin(theta) + yc;
xplot(:,1) = xT;
yplot(:,1) = yT;
for tstep = 1:ceil(tfinal/dt);
    for istep = 1:n
        a = xT(istep);
        b = yT(istep);
        ix = find(xv<=a,1,'last');
        jx = find(yv<=b,1,'last');
        uab = ...
            1/((xv(ix+1)-xv(ix))*(yv(jx+1)-yv(jx)))*...
            (u(ix,jx)*(xv(ix+1)-a)*(yv(jx+1)-b) + ...
            u(ix+1,jx)* (a-xv(ix)) * (yv(jx+1) - b) + ...
            u(ix,jx+1)*(xv(ix+1)-a)*(b-yv(jx)) + ...
            u(ix+1,jx+1)*(a-xv(ix))*(b-yv(jx)));
        vab = ...
            1/((xv(ix+1)-xv(ix))*(yv(jx+1)-yv(jx)))*...
            (v(ix,jx)*(xv(ix+1)-a)*(yv(jx+1)-b) + ...
            v(ix+1,jx)* (a-xv(ix)) * (yv(jx+1) - b) + ...
            v(ix,jx+1)*(xv(ix+1)-a)*(b-yv(jx)) + ...
            v(ix+1,jx+1)*(a-xv(ix))*(b-yv(jx)));
        xTn(istep,tstep) = xT(istep) + c*dt*uab;
        yTn(istep,tstep) = yT(istep) + c*dt*vab;
        oxplot(istep,tstep) = xv(ix);
        oyplot(istep,tstep) = yv(jx);
        quplot(istep,tstep) = vab;
        qvplot(istep,tstep) = uab;
    end
    xT = xTn(:,tstep)';
    yT = yTn(:,tstep)';
end

%% Plot of all time figure 1
q = 1:n/20:n;
xq = x(q,q);
yq = y(q,q);
uq =  u(q,q);
vq =  v(q,q);
figure(1);clf;
for ipp = 1:ceil(tfinal/dt);
    quiver(xq,yq,uq,vq)
    hold on;
    plot(xTn(:,ipp),yTn(:,ipp));
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    hold off;
    pause(0.001)
end
title('plot finished');
%% Plot of discrete times figure 2
n1 = 1;
n2 = round(totStep/4  );
n3 = round(totStep/2  );
n4 = round(totStep*3/4);
n5 = totStep;

t1 = dt*n1;
t2 = dt*n2;
t3 = dt*n3;
t4 = dt*n4;
t5 = dt*n5;
figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
plot(xTn(:,n1),yTn(:,n1),...
    xTn(:,n2),yTn(:,n2),...
    xTn(:,n3),yTn(:,n3),...
    xTn(:,n4),yTn(:,n4),...
    xTn(:,n5),yTn(:,n5))
hold on;
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
axis('square')
hLegend = legend(gca,...
    sprintf('t = %0.2f s',t1),...
    sprintf('t = %0.2f s',t2),...
    sprintf('t = %0.2f s',t3),...
    sprintf('t = %0.2f s',t4),...
    sprintf('t = %0.2f s',t5),'location','SouthEast');
hTitle = title(['Contour Plot of '...
    sprintf('$\\mathbf{\\phi}$ with velocity vectors')],...
    'Interpreter','latex');
set(gca, ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XMinorTick'  , 'on'          , ...
    'YMinorTick'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1             );
hold off;
printYesNo = 0;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('correctedLevelSetMethod')]);
    makeTable({'dt';'dx';'n'},[dt, dx, n]','myTable.tex')
end

