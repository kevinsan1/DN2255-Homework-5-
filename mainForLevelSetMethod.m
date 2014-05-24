clearvars;clc;close all force;
global n;
n = 50;
%% dt = 0.1 dx
printYesNo = 1;
dtp1 = 0.1;
dtplabel1 = '0.1\Deltax';
[dt,dx] = levelSetMethod(dtp1,dtplabel1);
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dx')]);
    makeTable({'dt';'dx';'n'},{dtplabel1;'L/(n-1)';'-'},[dt, dx, n]','levelSetTable1.tex')
end
%% dt = 1 dx
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
    makeTable({'dt';'dx';'n'},{dtplabel2;'L/(n-1)';'-'},[dt, dx, n]','levelSetTable2.tex')
end
%% dt = 1 dx
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
    makeTable({'dt';'dx';'n'},{'dx';'L/(n-1)';'-'},[dt, dx, n]','levelSetTable3.tex')
end
%% dt = 0.1 dx n = 200
clearvars -except n;clc;close all force;
n = 200;
dtp4 = 0.1 ;
dtplabel4 = '0.1\Deltax';
[dt,dx] = levelSetMethod(dtp4,dtplabel4);
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dxn200')]);
    makeTable({'dt';'dx';'n'},{'0.1dx';'L/(n-1)';'-'},[dt, dx, n]'...
        ,'levelSetTable4.tex')
end