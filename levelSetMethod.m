%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
tic
clear all;clc;close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 100;
dx = 2/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-1:dx:1,-1:dy:1);
% g = 2:n+1; % ghost cell index
%% Define velocity field
% u=-cos(pi*(X+1/2))*sin(3*pi/8*Y);
% v=sin(pi*(X+1/2))*cos(3*pi/8*Y);
u=2-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=2+sin(pi*(X+0.5)).*cos(3*pi/8*Y);
tFinal = 1.6;
dt = .08*dx;
tSteps = ceil(tFinal/dt);
%% Signed Distance Function
phi = signedDistance;
% phin = zeros(n+2,n+2);

%% Add Ghost Cells and Set Values
% g = 2:n+1 - ghost cell index
zeroBC = 'zeros';
rbc = 'reflective velocity BC';
pbc = 'periodic BC';
nfbc = 'no flux BC';
velBC = zeros;
phiBC = nfbc;
% u = setBoundariesOfNxN(u,velBC);
% v = setBoundariesOfNxN(v,velBC);
g=2:n+1;
up = zeros(n+2);
um = zeros(n+2);
vp = zeros(n+2);
vm = zeros(n+2);
up(g,g) = u(:,:);
um(g,g) = u(:,:);
vp(g,g) = v(:,:);
vm(g,g) = v(:,:);
%%
g = 2:n-1;
% phi = setBoundariesOfNxN(phi,phiBC);
%% Set u and v
% im = [n,[1:n-1]];
% ip = [[2:n],1];
% up = u(im,:);
% um = u(ip,:);
up(up<=0) = 0;
um(um>=0) = 0;
% vp = v(:,im);vm = v(:,ip);
vp(vp<=0) = 0;
vm(vm>=0) = 0;
myHandle = waitbar(0,'Initializing waitbar...');
tic;
%% Plotting variables
nplots = 100;        % Desired number of plots
plotStep = tSteps/nplots; % Number of steps between plots
iplot = 1;
tplot(1) = 0;
%% Main loop

for tn = 1:tSteps
    for i = 2:n-1
        for j = 2:n-1
            wxm = phi(i,j,tn) - phi(i-1,j,tn);	% x backward difference
            wxp = phi(i+1,j,tn) - phi(i,j,tn); 	% x forward difference
            wym = phi(i,j,tn) - phi(i,j-1,tn); 	% y backward difference
            wyp = phi(i,j+1,tn) - phi(i,j,tn); 	% y forward difference
            phi(i,j,tn+1) = phi(i,j,tn)...
                -dt/dx * ( up(i,j)*wxm + um(i,j)*wxp )...
                -dt/dy * ( vp(i,j)*wym + vm(i,j)*wyp );
        end
    end
    % phin = setBoundariesOfNxN(phin,phiBC);
    phi(1,:,tn+1)       = 	phi(n-1,:,tn+1);
    phi(n,:,tn+1)     = 	phi(2,:,tn+1);
    phi(:,1,tn+1)       = 	phi(:,n-1,tn+1);
    phi(:,n,tn+1)     =	phi(:,2,tn+1);
    if( rem(tSteps,plotStep) < 1 )  % Every plot_iter steps record
        iplot = iplot+1;
        phiplot(:,:,iplot) = phi(:,:,tn+1);       % Record a(i) for ploting
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
qp = 1:n/25:n;
for iPlot = 1:16:tSteps
    figure(2)
    contour(X,Y,phi(:,:,iPlot),[0,0],'r');
    hold on;
    contour(X,Y,phi(:,:,1),[0,0],'m');
    quiver(X(qp,qp),Y(qp,qp),u(qp,qp),v(qp,qp))
    hold off
    axis([-1 1 -1 1])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
    pause(.01);
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
% figure(4)
% contour(X,Y,phi(g,g,1),[0,0],'r'); % initial state
% hold on;
% contour(X,Y,phi(g,g,iplot),[0,0],'m'); % final state
% quiver(X(qp,qp),Y(qp,qp),u(g(qp),g(qp)),v(g(qp),g(qp)))
% hold off
% axis([-1 1 -1 1])
% axis('square')
% title(sprintf('t = %0.3g',iplot*dt));
%%
% savePath = ['/Users/kevin/SkyDrive/'...
%     'KTH Work/Period 3 2014/DN2255/Homework/5/matlab/mat/'];
% save([savePath sprintf('dt_1-0050e-4',dt)],...
%     'X','Y','u','v','x','y','tSteps','phi','dt')
% save([savePath sprintf('dt1-0050e-4phiplot',dt)],...
%     'X','Y','u','v','x','y','tSteps','phiplot','dt')
