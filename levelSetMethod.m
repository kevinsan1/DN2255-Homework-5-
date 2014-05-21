%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
tic
clearvars -except uother v other;clc;close all force;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
global n xc yc r;
%% Constants
n = 50;
xc = 0;
yc = -.6;
L = 2;
dx = L/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:dx:L/2,-L/2:dy:L/2);
tFinal = 1.6;
dt = .2*dx;
tSteps = ceil(tFinal/dt);
ia = 1:n;
im = [n,1:n-1];
ip = [2:n,1];
%% Signed Distance Function
phi = signedDistance;
%% Define velocity field
%     u=.8+cos(2*pi*Y);
%     v=.8+sin(2*pi*X);
u=2-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=2+sin(pi*(X+0.5)).*cos(3*pi/8*Y);
up = u;
um = u;
vp = v;
vm = v;
up(up<=0) = 0;
um(um>=0) = 0;
vp(vp<=0) = 0;
vm(vm>=0) = 0;

%% Plotting variables
nplots = 20;        % Desired number of plots
plotStep = round(tSteps/nplots); % Number of steps between plots
iplot = 1;
tplot(1) = 0;
phiplot(:,:,1) = phi;
%% Initialize variables
fx = zeros(n);
fy = zeros(n);
phinew = zeros(n);
tplot = zeros(1,nplots+1);
phiplot = zeros(n,n,nplots+1);
%% Main loop
myHandle = waitbar(0,'Initializing waitbar...');
tic;

for tn = 1:tSteps
    for i = 1:n
        for j = 1:n
            wxm = phi(i,j) - phi(im(i),j);	% x backward difference
            wxp = phi(ip(i),j) - phi(i,j); 	% x forward difference
            wym = phi(i,j) - phi(i,im(j)); 	% y backward difference
            wyp = phi(i,ip(j)) - phi(i,j); 	% y forward difference
            fx(i,j) = up(im(i),j)*wxm + um(ip(i),j)*wxp;
            fy(i,j) = vp(i,im(j))*wym + vm(i,ip(j))*wyp;
            phinew(i,j) = phi(i,j)-dt*( fx(i,j)/dx + fy(i,j)/dy );
        end
    end
    % phin = setBoundariesOfNxN(phin,phiBC);
    
    phi = phinew;
    if( rem(tn,plotStep) < 1 )  % Every plot_iter steps record
        iplot = iplot+1;
        phiplot(:,:,iplot) = phi;       % Record a(i) for ploting
        tplot(iplot) = dt*tn;
    end
    %% Here's the progress bar code
    time=toc;
    Perc=tn/tSteps;
    Trem=time/Perc-time; %Calculate the time remaining
    Hrs=floor(Trem/3600);Min=floor((Trem-Hrs*3600)/60);
    waitbar(Perc,myHandle,[sprintf('%0.1f',Perc*100) '%, '...
        sprintf('%03.0f',Hrs) ':'...
        sprintf('%02.0f',Min) ':'...
        sprintf('%02.0f',rem(Trem,60)) ' remaining']);
end
%% Plotting
% [x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
% u=-cos(pi*(x+1/2))*sin(3*pi/8*y);
% v=sin(pi*(x+1/2))*cos(3*pi/8*y);
% u = 2+cos(2*pi*x).*sin(2*pi*y);
% v = 2+sin(2*pi*x);
%%
qp = 1:round(n/25):n;
% u=2-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
% v=2+sin(pi*(X+0.5)).*cos(3*pi/8*Y);
quivU = up(im(:),:) + um(ip(:),:);
quivV = vp(:,im(:)) + vm(:,ip(:));
figure(2);clf;
for iPlot = 1:iplot
    figure(2)
    a1=contour(X,Y,phiplot(:,:,iPlot),[0,0],'r');
    hold on;
    a2=contour(X,Y,phiplot(:,:,1),[0,0],'b');
    a3=quiver(X(qp,qp),Y(qp,qp),quivU(qp,qp),quivV(qp,qp));
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
    hold off
    pause(.001);
end
% %%

% for iPlot = 1:iplot
%     figure(2)
%     contour(X,Y,phi(:,:,iPlot),[0,0],'r');
%     hold on;
%     contour(X,Y,phi(:,:,1),[0,0],'m');
%     quiver(X(qp,qp),Y(qp,qp),1./u(g(qp),g(qp)),1./v(g(qp),g(qp)))
%     hold off
%     axis([-1 1 -1 1])
%     axis('square')
%     title(sprintf('t = %0.3g',iPlot*dt));
% end
%%
% for iPlot = 1:iplot
%     figure(3)
%     mesh(X,Y,phi(:,:,iPlot));
%     title(sprintf('t = %0.3g',iPlot*dt));
%     axis([-1 1 -1 1 min(min(phi(:,:,1))) 0]);
%     %     pause(0.1)
% end
%% Plot for print
% clc;clf;
% 
% quivU = up(im(:),:) + um(ip(:),:);
% quivV = vp(:,im(:)) + vm(:,ip(:));
% qp = 1:round(n/25):n;
% n1 = 1;
% n2 = 4;
% n3 = 6;
% n4 = 8;
% n5 = 9;
% 
% t1 = tplot(n1);
% t2 = tplot(n2);
% t3 = tplot(n3);
% t4 = tplot(n4);
% t5 = tplot(n5);
% figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);
% hold on;
% y1_plot = contour(X,Y,phiplot(:,:,n1),[0,0],'b'); % initial state
% y2_plot = contour(X,Y,phiplot(:,:,n2),[0,0],'m'); % second state
% y3_plot = contour(X,Y,phiplot(:,:,n3),[0,0],'g'); % third state
% y4_plot = contour(X,Y,phiplot(:,:,n4),[0,0],'k'); % fourth state
% y5_plot = contour(X,Y,phiplot(:,:,n5),[0,0],'r'); % final state
% quivP = quiver(X(qp,qp),Y(qp,qp),quivU(qp,qp),quivV(qp,qp),'b');
% axis([-L/2 L/2 -L/2 L/2])
% axis('square')
% hLegend = legend(gca,...
%     sprintf('t = %g s',t1),...
%     sprintf('t = %g s',t2),...
%     sprintf('t = %g s',t3),...
%     sprintf('t = %g s',t4),...
%     sprintf('t = %g s',t5),'location','Best');
% hTitle = title(sprintf('Contour Plot of $\\mathbf{\\phi}$ with velocity vectors'),'Interpreter','latex');
% hold off;
% set(gca, ...
%     'Box'         , 'off'         , ...
%     'TickDir'     , 'out'         , ...
%     'TickLength'  , [.02 .02]     , ...
%     'XMinorTick'  , 'on'          , ...
%     'YMinorTick'  , 'on'          , ...
%     'XColor'      , [.3 .3 .3]    , ...
%     'YColor'      , [.3 .3 .3]    , ...
%     'LineWidth'   , 1             );
% printYesNo = 0;
% if printYesNo == 1
%     saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
%         ' 3 2014/DN2255/Homework/5/Figures/'];
%     addpath(saveFigurePath);
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-depsc2', [saveFigurePath ...
%         sprintf('contourPlotOfQ1PartBbigtimestep')]);
%     makeTable({'dt';'dx';'n'},[dt, dx, n]','myTable.tex')
% end
%%
% savePath = ['/Users/kevin/SkyDrive/'...
%     'KTH Work/Period 3 2014/DN2255/Homework/5/matlab/mat/'];
% save([savePath sprintf('dt_1-0050e-4',dt)],...
%     'X','Y','u','v','x','y','tSteps','phi','dt')
% save([savePath sprintf('dt1-0050e-4phiplot',dt)],...
%     'X','Y','u','v','x','y','tSteps','phiplot','dt')
