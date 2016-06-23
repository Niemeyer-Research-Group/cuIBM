close all
clear
clc

M = dlmread('/scratch/src/cuIBM/validation/luo/test/forces','\t');

y = fft(M(40:end,2))
plot(M(40:end,1),y)