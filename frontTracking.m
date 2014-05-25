clearvars;clc;
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
dt = 1*dx;
tfinal = 2;
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

for tstep = 1:totStep;
    for istep = 2:n-1
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
    xTn(1,tstep) = xTn(n-1,tstep);
    yTn(1,tstep) = yTn(n-1,tstep);
    xTn(n,tstep) = xTn(2,tstep);
    yTn(n,tstep) = yTn(2,tstep);
    xT = xTn(:,tstep)';
    yT = yTn(:,tstep)';
%     quiver(xq,yq,uq,vq)
%     hold on;
%     plot(xTn(:,tstep),yTn(:,tstep));
%     axis([-L/2 L/2 -L/2 L/2])
%     axis('square')
%     hold off
end

%% Plot of all time figure 1
q = 1:n/20:n;
xq = x(q,q);
yq = y(q,q);
uq =  u(q,q);
vq =  v(q,q);
figure(1);clf;
for ipp = 1:4:tstep-1;
    quiver(xq,yq,uq,vq)
    hold on;
    plot(xTn(:,ipp),yTn(:,ipp));
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    hold off;
    pause(0.1);
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
figure1;
hold on;
plot( ...
	xTn(:,n1),yTn(:,n1),...
    xTn(:,n2),yTn(:,n2),...
    xTn(:,n3),yTn(:,n3),...
    xTn(:,n4),yTn(:,n4),...
    xTn(:,n5),yTn(:,n5))
hold on;
hXLabel = xlabel('x');
hYLabel = ylabel('y');
quiver(xq,yq,uq,vq)
axis([-L/2 L/2 -L/2 L/2])
axis('square')
hText = text(-0.95,0,...
    [sprintf('N = %g\n\\Deltat = 0.1\\Deltax',n) sprintf(' = %4.4f\n\\Deltax = L/(N-1) = %4.4f',dt,dx)],...
    'EdgeColor',[0 0 0],...
    'BackgroundColor',[1 1 1]);
hLegend = legend(gca,...
    sprintf('t = %0.2f s',t1),...
    sprintf('t = %0.2f s',t2),...
    sprintf('t = %0.2f s',t3),...
    sprintf('t = %0.2f s',t4),...
    sprintf('t = %0.2f s',t5),'location','SouthEast');
hTitle = title(['Contour Plot of '...
    sprintf('$\\mathbf{\\phi}$ with velocity vectors')],...
    'Interpreter','latex');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel, hText], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel, hText]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );
set(gca, ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XMinorTick'  , 'on'          , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on',...
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
        sprintf('frontTracking')]);
end

