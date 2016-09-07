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

%% plot force of osc cylinders with no flow
% close all
figure
clear
clc
path = '/scratch/src/cuIBM/validation/osc/static/forces';
delim = '\t';

%read data
forces = dlmread(path,delim,1,0);

%plot
hold on
plot(forces(10:end,1),forces(10:end,2),'r')
hold off
xlabel('time')
ylabel('force')
title('Oscillating Cylinder in Static Flow')
axis([0 2 -0.5 0.5])

%% plot force of osc cylinders with flow
% close all
figure
clear
clc
path = '/scratch/src/cuIBM/validation/osc/flow/forces';
delim = '\t';

%read data
forces = dlmread(path,delim,1,0);
%plot
hold on
plot(forces(10:end,1),forces(10:end,2),'r')
hold off
xlabel('time')
ylabel('force')
title('Drag')
axis([0 10 0 2])

%% plot foce of all osc cylinders with flow

% close all
% clear
% clc

pathbase = '/scratch/src/cuIBM/validation/osc/flow/';
delim = '\t';

%read data
a = dlmread(strcat(pathbase,'a','/forces'),delim,1,0);
b = dlmread(strcat(pathbase,'b','/forces'),delim,1,0);
c = dlmread(strcat(pathbase,'c','/forces'),delim,1,0);
d = dlmread(strcat(pathbase,'d','/forces'),delim,1,0);
e = dlmread(strcat(pathbase,'e','/forces'),delim,1,0);
f = dlmread(strcat(pathbase,'f','/forces'),delim,1,0);
g = dlmread(strcat(pathbase,'g','/forces'),delim,1,0);
h = dlmread(strcat(pathbase,'h','/forces'),delim,1,0);
i = dlmread(strcat(pathbase,'i','/forces'),delim,1,0);
j = dlmread(strcat(pathbase,'j','/forces'),delim,1,0);
k = dlmread(strcat(pathbase,'k','/forces'),delim,1,0);
l = dlmread(strcat(pathbase,'l','/forces'),delim,1,0);

hold on
%row 1
subplot(4,3,1)
plot(a(:,1),a(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,3,2)
plot(b(:,1),b(:,2),'k')
axis([0 10 0 2])

subplot(4,3,3)
plot(c(:,1),c(:,2),'k')
axis([0 10 0 2])

%row 2
subplot(4,3,4)
plot(d(:,1),d(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,3,5)
plot(e(:,1),e(:,2),'k')
axis([0 10 0 2])

subplot(4,3,6)
plot(f(:,1),f(:,2),'k')
axis([0 10 0 2])

%row 3
subplot(4,3,7)
plot(g(:,1),g(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,3,8)
plot(h(:,1),h(:,2),'k')
axis([0 10 0 2])

subplot(4,3,9)
plot(i(:,1),i(:,2),'k')
axis([0 10 0 2])

%row 4
subplot(4,3,10)
plot(j(:,1),j(:,2),'k')
ylabel('force')
xlabel('time')
axis([0 10 0 2])

subplot(4,3,11)
plot(k(:,1),k(:,2),'k')
xlabel('time')
axis([0 10 0 2])

subplot(4,3,12)
plot(l(:,1),l(:,2),'k')
xlabel('time')
axis([0 10 0 2])
hold off
