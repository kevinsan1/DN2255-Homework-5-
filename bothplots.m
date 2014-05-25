clearvars -except n;clc;
global n;
n = 200;
dtp4 = 0.1 ;
dtplabel4 = '0.1\Deltax';
[dt,dx,phi,X,Y,tplot,figure1] = levelSetMethod(dtp4,dtplabel4);
%% Plot for print level set
u=-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=sin(pi*(X+0.5)).*cos(3*pi/8*Y);
L=2;
qp = 1:round(n/20):n;
totTime=length(tplot);
n1 = 1;
n2 = round(1/4*totTime);
n3 = round(.5*totTime);
n4 = round(3/4*totTime);
n5 = totTime;

t1 = tplot(n1);
t2 = tplot(n2);
t3 = tplot(n3);
t4 = tplot(n4);
t5 = tplot(n5);
figure1=figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);hold on;
[C1,h] = contour(X,Y,phi(:,:,n1),[0,0],'r'); % initial state
[C2,h] = contour(X,Y,phi(:,:,n2),[0,0],'r'); % second state
[C3,h] = contour(X,Y,phi(:,:,n3),[0,0],'r'); % third state
[C4,h] = contour(X,Y,phi(:,:,n4),[0,0],'r'); % fourth state
[C5,h] = contour(X,Y,phi(:,:,n5),[0,0],'r'); % final state
hold off;
%%
n = 200;
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
myHandle = waitbar(0,'Initializing waitbar...');
tic;
for tstep = 1:totStep;
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
        %% Here's the progress bar code
    time=toc;
    Perc=tstep/totStep;
    Trem=time/Perc-time; %Calculate the time remaining
    Hrs=floor(Trem/3600);Min=floor((Trem-Hrs*3600)/60);
    waitbar(Perc,myHandle,[sprintf('%0.1f',Perc*100) '%, '...
        sprintf('%03.0f',Hrs) ':'...
        sprintf('%02.0f',Min) ':'...
        sprintf('%02.0f',rem(Trem,60)) ' remaining front method']);
end
%%
x1 = C1(1,2:end);
y1 = C1(2,2:end);
x2 = C2(1,2:end);
y2 = C2(2,2:end);
x3 = C3(1,2:end);
y3 = C3(2,2:end);
x4 = C4(1,2:end);
y4 = C4(2,2:end);
x5 = C5(1,2:end);
y5 = C5(2,2:end);
figure1=figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);hold on;
plot(xTn(:,n1),yTn(:,n1),'b',x1,y1,'r');
plot( ...
    xTn(:,n2),yTn(:,n2),'b',...
    xTn(:,n3),yTn(:,n3),'b',...
    xTn(:,n4),yTn(:,n4),'b',...
    xTn(:,n5),yTn(:,n5),'b');
plot(x2,y2,'r',x3,y3,'r',x4,y4,'r',x5,y5,'r');
quivP = quiver(X(qp,qp),Y(qp,qp),u(qp,qp),v(qp,qp),...
    'color',[.5 .5 .5]);
axis([-L/2 L/2 -L/2 L/2])
axis('square')
hXLabel = xlabel('x');
hYLabel = ylabel('y');
% Create axes
hText = text(-0.95,0,...
    [sprintf('N = %g\n\\Deltat = 0.1',n) sprintf(' = %4.4f\n\\Deltax = L/(N-1) = %4.4f',dt,dx)],...
    'EdgeColor',[0 0 0],...
    'BackgroundColor',[1 1 1]);
%
hLegend = legend(...
    sprintf('Front Tracking'),...
    sprintf('Level Set'),...
    'location','SouthEast');
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
    'YGrid'     ,'on',...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1             );
hold off;
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('frontandLevelSet')]);
end