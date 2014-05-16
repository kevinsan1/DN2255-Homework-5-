% clear all;clc;close all;
clearvars -except phi
a = magic(5)
b = zeros(7)
b(2:6,2:6) = a;
b(1,:) = b(2,:)
b(:,1) = b(:,2)
b(end,:) = b(end-1,:)
b(:,end)=b(:,end-1)
