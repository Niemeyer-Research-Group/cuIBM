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

%% plot force of osc cylinders with flow
close all
clear
clc
pathbase = '/scratch/src/cuIBM/validation/';
cases = 'osc/flow/forces';
delim = '\t';

%read data
forces = dlmread(strcat(pathbase,cases),delim,0,0);
x = dlmread(strcat(pathbase,'osc/flow/midPosition'),delim,1,0);
%plot
hold on
plot(forces(10:end,1),forces(10:end,2),'r') %present
hold off
xlabel('time')
ylabel('force')
title('Drag')
axis([0 10 -2 6])

% figure
% X = 0:0.1:10;
% plot(forces(10:end,1),x(10:end,2),'k'), hold on
% plot(X,-0.25*cos(2*pi*0.2*X),'ko')
% plot(forces(10:end,1),x(10:end,4),'r')
% plot(X,0.1*pi*sin(2*pi*0.2*X),'ro')
% xlabel('time')
% legend('position','expected position','velocity','expected velocity')
% ylabel('x position')

%% plot foce of all osc cylinders with flow

% close all
% clear
% clc

pathbase = '/scratch/src/cuIBM/validation/osc/';
cases = '/flow/forces';
delim = '\t';

%read data
a = dlmread(strcat(pathbase,'ab','/forcesa'),delim,0,0);
b = dlmread(strcat(pathbase,'ab','/forcesb'),delim,0,0);
c = dlmread(strcat(pathbase,'cd','/forcesc'),delim,0,0);
d = dlmread(strcat(pathbase,'cd','/forcesd'),delim,0,0);
e = dlmread(strcat(pathbase,'ef','/forcese'),delim,0,0);
f = dlmread(strcat(pathbase,'ef','/forcesf'),delim,0,0);
g = dlmread(strcat(pathbase,'gh','/forcesg'),delim,0,0);
h = dlmread(strcat(pathbase,'gh','/forcesh'),delim,0,0);

%plot
hold on
% subplot(4,2,1)
% plot(a(:,1),a(:,2),'k')
% ylabel('force')
% axis([0 10 -2 6])
% 
% subplot(4,2,2)
% plot(b(:,1),b(:,2),'k')
% axis([0 10 -2 6])
% 
% subplot(4,2,3)
% plot(c(:,1),c(:,2),'k')
% ylabel('force')
% axis([0 10 -2 6])
% 
% subplot(4,2,4)
% plot(d(:,1),d(:,2),'k')
% axis([0 10 -2 6])
% 
% subplot(4,2,5)
% plot(e(:,1),e(:,2),'k')
% ylabel('force')
% axis([0 10 -2 6])
% 
% subplot(4,2,6)
% plot(f(:,1),f(:,2),'k')
% axis([0 10 -2 6])

% subplot(4,2,7)
plot(g(:,1),g(:,2),'k')
% xlabel('time')
% ylabel('force')
% axis([0 10 -2 6])

% subplot(4,2,8)
% plot(h(:,1),h(:,2),'r')
% xlabel('time')
% axis([0 10 -2 6])

hold off
