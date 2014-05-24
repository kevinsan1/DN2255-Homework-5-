clearvars;clc;close all force;
global n;
n = 50;
%% dt = 0.1 dx
printYesNo = 1;
dtp1 = 0.1;
dtplabel1 = '0.1dx';
[dt,dx] = levelSetMethod(dtp1,dtplabel1);
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt0p1dx')]);
    makeTable({'dt';'dx';'n'},{dtplabel1;'L/(n-1)';'-'},[dt, dx, n]','levelSetTable.tex')
end
%% dt = 1 dx
clearvars -except n;clc;close all force;
dtp2 = 1 ;
dtplabel2 = 'dx';
[dt,dx] = levelSetMethod(dtp2,dtplabel2);
printYesNo = 1;
if printYesNo == 1
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/Period' ...
        ' 3 2014/DN2255/Homework/5/Figures/'];
    addpath(saveFigurePath);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('dt1dx')]);
    makeTable({'dt';'dx';'n'},{dtplabel2;'L/(n-1)';'-'},[dt, dx, n]','levelSetTable.tex')
end