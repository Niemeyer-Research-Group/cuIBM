%plot pressure
clc
clear all
close all
number = '100';
type = 'p';
path = strcat('/scratch/src/cuIBM/validation/luo/test/output/',number,type,'.csv');
delim = '\t';
M = dlmread(path,delim,1,0);
surf(M)
xlabel('x')
ylabel('y')