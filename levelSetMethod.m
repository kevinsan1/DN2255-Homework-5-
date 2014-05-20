%% Level Set Method
% phi_t + u_{ext} * phi_x = 0
tic
clearvars -except uother v other;clc;close all force;
myPath = ['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 3 2014/DN2255/Homework/5/matlab/'];
cd(myPath);
addpath(genpath(myPath));
global n; global r; global insideCircleTest;
%% Constants
n = 100;
L = 2;
dx = L/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:dx:L/2);
tFinal = 1.6;
dt = .1*dx;
tSteps = ceil(tFinal/dt);
ia = 1:n;
im = [n,1:n-1];
ip = [2:n,1];
%% Signed Distance Function
phi = signedDistance;
%% Define velocity field
%     u=.8+cos(2*pi*Y);
%     v=.8+sin(2*pi*X);
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
plotStep = tSteps/nplots; % Number of steps between plots
iplot = 1;
tplot(1) = 0;
phiplot(:,:,1) = phi;
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
            fx(i,j) = max(u(im(i),j),0)*wxm + min(u(ip(i),j),0)*wxp;
            fy(i,j) = max(v(i,im(j)),0)*wym + min(v(i,ip(j)),0)*wyp;
            phinew(i,j) = phi(i,j)-dt*( fx(i,j)/dx + fy(i,j)/dy );
        end
    end
    % phin = setBoundariesOfNxN(phin,phiBC);
    phi = phinew;
    if( rem(tSteps,plotStep) < 1 )  % Every plot_iter steps record
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
for iPlot = 1:iplot
    figure(2)
    contour(X,Y,phiplot(:,:,iPlot),[0,0],'r');
    hold on;
    contour(X,Y,phiplot(:,:,1),[0,0],'b');
    quiver(X(qp,qp),Y(qp,qp),quivU(qp,qp),quivV(qp,qp))
    hold off
    axis([-L/2 L/2 -L/2 L/2])
    axis('square')
    title(sprintf('t = %0.3g',iPlot*dt));
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
