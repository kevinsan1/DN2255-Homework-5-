function [dt,dx] = levelSetMethod(dtfraction,myTextLabel)
%  	Description
%	out = levelSetMethod(n,dt,none,none)
%   n = number of grids
%	dt = 
%	none =
%	none =
global n xc yc r;
% Long description
if nargin < 1
	L = 2;
	dx = L/(n-1);
	dt = 0.1*dx;
end
%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
%% Constants
xc = 0;
yc = -.6;
r = 0.3;
L = 2;
dx = L/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:dx:L/2,-L/2:dy:L/2);
tFinal = 1.6;
dt = dtfraction*dx;
tSteps = ceil(tFinal/dt);
ia = 1:n;
im = [n,1:n-1];
ip = [2:n,1];
%% Signed Distance Function
phi = signedDistance;
%% Define velocity field
u=-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=sin(pi*(X+0.5)).*cos(3*pi/8*Y);
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
phiplot = phi;
%% Main loop
for tn = 1:tSteps
    for i = 2:n-1
        for j = 2:n-1
            wxm = phi(i,j) - phi(i-1,j);	% x backward difference
            wxp = phi(i+1,j) - phi(i,j); 	% x forward difference
            wym = phi(i,j) - phi(i,j-1); 	% y backward difference
            wyp = phi(i,j+1) - phi(i,j); 	% y forward difference
            fx(i,j) = vp(i,j)*wxm + vm(i,j)*wxp;
            fy(i,j) = up(i,j)*wym + um(i,j)*wyp;
            phinew(i,j) = phi(i,j)-dt*( fx(i,j)/dx + fy(i,j)/dy );
        end
    end
    phi = phinew;
    if( rem(tn,plotStep) < 1 ) % Every plot_iter steps record
        iplot = iplot+1;
        phiplot(:,:,iplot) = phi; % Record phi(i) for ploting
        tplot(iplot) = dt*tn;
    end
end

%%
% qp = 1:round(n/25):n;
% % u=2-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
% % v=2+sin(pi*(X+0.5)).*cos(3*pi/8*Y);
% quivU = up(im(:),:) + um(ip(:),:);
% quivV = vp(:,im(:)) + vm(:,ip(:));
% figure(1);clf;
% for iPlot = 1:iplot
%     figure(1)
%     a1=contour(X,Y,phiplot(:,:,iPlot),[0,0],'r');
%     hold on;
%     a2=contour(X,Y,phiplot(:,:,1),[0,0],'b');
%     a3=quiver(X(qp,qp),Y(qp,qp),quivU(qp,qp),quivV(qp,qp));
%     axis([-L/2 L/2 -L/2 L/2])
%     axis('square')
%     title(sprintf('t = %0.3g',iPlot*dt));
%     hold off
%     pause(.001);
% end

%% Plot for print

quivU = up(im(:),:) + um(ip(:),:);
quivV = vp(:,im(:)) + vm(:,ip(:));
qp = 1:round(n/25):n;
n1 = 1;
n2 = 5;
n3 = 10;
n4 = 15;
n5 = 20;

t1 = tplot(n1);
t2 = tplot(n2);
t3 = tplot(n3);
t4 = tplot(n4);
t5 = tplot(n5);
figure1=figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);clf;hold on;
quivP = quiver(X(qp,qp),Y(qp,qp),quivU(qp,qp),quivV(qp,qp),...
    'color',[.5 .5 .5]);
y1_plot = contour(X,Y,phiplot(:,:,n1),[0,0],'r'); % initial state
y2_plot = contour(X,Y,phiplot(:,:,n2),[0,0],'m'); % second state
y3_plot = contour(X,Y,phiplot(:,:,n3),[0,0],'color',[.2 .5 .9]); % third state
y4_plot = contour(X,Y,phiplot(:,:,n4),[0,0],'b'); % fourth state
y5_plot = contour(X,Y,phiplot(:,:,n5),[0,0],'color',[.3 .3 .3]); % final state
axis([-L/2 L/2 -L/2 L/2])
axis('square')
hXLabel = xlabel('x');
hYLabel = ylabel('y');
% Create axes
hText = text(-0.95,0,...
    sprintf(['n = %g\ndt = ' myTextLabel ' = %4.4f\ndx = L/(N-1) = %4.4f'],n,dt,dx),...
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
  'YGrid'       , 'on'          , ...
  'XColor'      , [.3 .3 .3]    , ...
  'YColor'      , [.3 .3 .3]    , ...
  'LineWidth'   , 1             );
hold off;
%%
% savePath = ['/Users/kevin/SkyDrive/'...
%     'KTH Work/Period 3 2014/DN2255/Homework/5/matlab/mat/'];
% save([savePath sprintf('dt_1-0050e-4',dt)],...
%     'X','Y','u','v','x','y','tSteps','phi','dt')
% save([savePath sprintf('dt1-0050e-4phiplot',dt)],...
%     'X','Y','u','v','x','y','tSteps','phiplot','dt')
end % function
