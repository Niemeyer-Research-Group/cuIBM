clear
close
clc

%% clyiner
%setup filepaths
pathbase = '/scratch/src/cuIBM/validation/';
cases = 'cylinder/Re40/forces';
validationPath = '/scratch/src/cuIBM/validation-data/cylinderRe40-KL95.txt';
delim = '\t';

%read data
forces = dlmread(strcat(pathbase,cases),delim,1,0);
validation = dlmread(validationPath,delim,0,0);

%plot
hold on
plot(forces(:,1),forces(:,2)*2,'r') %present
plot(validation(:,1)*0.5,validation(:,2),'ko') %validation
% plot(forces(20:end,1),forces(20:end,3),'r') %pressure
% plot(forces(20:end,1),forces(20:end,4),'b') %dudn
% plot(forces(20:end,1),forces(20:end,5)*2,'k') %sum
% plot(forces(20:end,1),forces(20:end,2)*2,'r') %og
% plot(forces(20:end,1),forces(20:end,5)*2,'b') %mine
% plot(validation(:,1)*0.5,validation(:,2),'ko') %validation

hold off
%beautify
legend('Present Work','Validation')
% legend('pre','dudn','sum')
% legend('og','new','validation')
xlabel('time')
ylabel('force')
title('Drag')
axis([0 5 0 6])

%% 
close all
pathbase = '/scratch/src/cuIBM/validation/';
cases = 'osc/flow/forces';
delim = '\t';

%read data
forces = dlmread(strcat(pathbase,cases),delim,0,0);
x = dlmread(strcat(pathbase,'osc/flow/midPosition'),delim,1,0);
%plot
hold on
plot(forces(10:end,1)*5,forces(10:end,2)*5,'r') %present
hold off
xlabel('time')
ylabel('force')
title('Drag')
axis([0 10 -2 6])

figure
plot(forces(10:end,1)*5,x(10:end,2)), hold on
plot(forces(10:end,1)*5,x(10:end,4))
xlabel('time')
legend('position','velocity')
% ylabel('x position')

