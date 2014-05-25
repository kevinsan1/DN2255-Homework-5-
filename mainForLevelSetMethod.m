clearvars;clc;close all force;
addpath(genpath(['/Users/kevin/SkyDrive/KTH Work/' ...
'Period 3 2014/DN2255/Homework/5/matlab/functions']));
cd(['/Users/kevin/SkyDrive/KTH Work/Period' ...
' 3 2014/DN2255/Homework/5/matlab/']);
global n;
n = 50;
%% dt = 0.1 dx n = 50
printYesNo = 1;
dtp1 = 0.1;
dtplabel1 = '0.1\Deltax';
[dt,dx,phi,X,Y] = levelSetMethod(dtp1,dtplabel1);
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dx')]);
    makeTable({'dt';'dx';'n'},{dtplabel1;'L/(n-1)';'-'}, ...
    [dt, dx, n]','levelSetTable1.tex')
end
save t1p1
%% dt = 1 dx n = 50
clearvars -except n;clc;close all force;
dtp2 = 1 ;
dtplabel2 = '\Deltax';
[dt,dx] = levelSetMethod(dtp2,dtplabel2);
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt1dx')]);
    makeTable({'dt';'dx';'n'},{dtplabel2;'L/(n-1)';'-'}, ...
    [dt, dx, n]','levelSetTable2.tex')
end
%% dt = 1 dx n = 200
clearvars -except n;clc;close all force;
n = 200;
dtp3 = 1 ;
dtplabel3 = '\Deltax';
[dt,dx] = levelSetMethod(dtp3,dtplabel3);

printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt1dxn200')]);
    makeTable({'dt';'dx';'n'},{'dx';'L/(n-1)';'-'}, ...
    [dt, dx, n]','levelSetTable3.tex')
end
%% dt = 0.1 dx n = 200
clearvars -except n;clc;
n = 200;
dtp4 = 0.1 ;
dtplabel4 = '0.1\Deltax';
[dt,dx,phi,X,Y,tplot] = levelSetMethod(dtp4,dtplabel4);
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dxn200')]);
    makeTable({'dt';'dx';'n'},{'0.1dx';'L/(n-1)';'-'}, ...
    [dt, dx, n]','levelSetTable4.tex')
end
%% Check if a signed distance function
clf;
L = 2;
[C,h]=contour(X,Y,phi(:,:,end),...
    [0,-0.4],'ShowText','on')
hCLabel=clabel(C,h,'FontSize',10,'Color','k','Rotation',0,'LabelSpacing',100)
axis([-L/2 L/2 -L/2 L/2])
%
axis('square')
hXLabel = xlabel('x');
hYLabel = ylabel('y');
% Create axes
hTitle = title(['Contour Plot of '...
    sprintf('$\\mathbf{\\phi}$ with 2 contour levels')],...
    'Interpreter','latex');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set(gca             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
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
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1             );
hold off;

%% dt = 0.1 dx n = 200 and forward method
clearvars -except n;clc;
n = 200;
dtp4 = 0.1 ;
dtplabel4 = '0.1\Deltax';
[dt,dx,phi,X,Y,tplot,figure1] = levelSetMethod(dtp4,dtplabel4);
%%
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dxn200')]);
    makeTable({'dt';'dx';'n'},{'0.1dx';'L/(n-1)';'-'}, ...
    [dt, dx, n]','levelSetTable4.tex')
end
