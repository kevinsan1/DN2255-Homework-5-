function [mat1,mat2,mat3,mat4,mat5,mat6,mat7] = test(plotmat)
%script here
close all;clc;clearvars -except plotmat
myPath = ['/Users/kevin/SkyDrive',...
    '/KTH Work/Period 3 2014/DN2255/Homework/5/matlab/mat/'];
files = dir(fullfile(myPath, '*.mat'));
%%
[mat1,mat2,mat3,mat4,mat5,mat6,mat7] = files(:,1).name;
%%
if nargin<1
    plotmat = mat3;
end

f(plotmat);
end
%%
function  f(matfile)
%function here
load(matfile);
for iPlot = 1:15921
    figure(2)
    contour(X,Y,phiplot(2:end-1,2:end-1,iPlot),[0,0],'r');
    hold on;
    contour(X,Y,phiplot(2:end-1,2:end-1,1),[0,0],'g');
    quiver(x,y,u,v);
    hold off;
    axis([-1 1 -1 1]);
    axis('square');
    title([matfile sprintf('     t = %0.3g',iPlot*dt)]);
    pause(.0005);
end
end