%% Level Set Method
clear all;close all force;clc
%% Constants
global n xc yc r insideCircleTest;
n = 100;
xc = 0;
yc = -.6;
r = 0.3;
L = 2;
dx = L/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:dx:L/2,-L/2:dy:L/2);
tFinal = .5;
dt = .01*dx;
tSteps = ceil(tFinal/dt);
%% Boundary Conditions
im = [n,1:n-1]; % Periodic i-1
ip = [2:n,1]; % Periodic i+1
%% Signed Distance Function
phi = signedDistance; % signedDistance
%% Define velocity field
u=-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=sin(pi*(X+0.5)).*cos(3*pi/8*Y);
up = u;
um = u;
vp = v;
vm = v;
up(up<=0) = 0; % if up(index)<=0, set up(index) = 0
um(um>=0) = 0; % if um(index)>=0, set um(index) = 0
vp(vp<=0) = 0; % if vp(index)<=0, set vp(index) = 0
vm(vm>=0) = 0; % if um(index)>=0, set vm(index) = 0
%% Initialize variables
nplots = 20;        % number of plots
plotStep = round(tSteps/nplots); % Number of steps between plots
iplot = 1;
tplot(1) = 0;
phiplot = zeros(n,n,nplots+1);
phiplot(:,:,1) = phi;
fx = zeros(n);
fy = zeros(n);
phinew = zeros(n);
tplot = zeros(1,nplots+1);
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
    phi = phinew;
    if( rem(tn,plotStep) < 1 )    % Every plot iter steps record
        iplot = iplot+1;
        phiplot(:,:,iplot) = phi; % Record phi for ploting
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
%%
qp = 1:round(n/25):n;
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