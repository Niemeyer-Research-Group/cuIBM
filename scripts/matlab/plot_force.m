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
legend('Fx','FxX','FxY','FxU')
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
m = dlmread(strcat(pathbase,'m','/forces'),delim,1,0);
n = dlmread(strcat(pathbase,'n','/forces'),delim,1,0);
o = dlmread(strcat(pathbase,'o','/forces'),delim,1,0);
p = dlmread(strcat(pathbase,'p','/forces'),delim,1,0);
q = dlmread(strcat(pathbase,'q','/forces'),delim,1,0);
r = dlmread(strcat(pathbase,'r','/forces'),delim,1,0);
s = dlmread(strcat(pathbase,'s','/forces'),delim,1,0);
t = dlmread(strcat(pathbase,'t','/forces'),delim,1,0);

hold on
%row 1
subplot(4,2,1)
plot(a(:,1),a(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,2,2)
plot(e(:,1),e(:,2),'k')
axis([0 10 0 2])

%row 2
subplot(4,2,3)
plot(f(:,1),f(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,2,4)
plot(j(:,1),j(:,2),'k')
axis([0 10 0 2])

%row 3
subplot(4,2,5)
plot(k(:,1),k(:,2),'k')
ylabel('force')
axis([0 10 0 2])

subplot(4,2,6)
plot(o(:,1),o(:,2),'k')
axis([0 10 0 2])

%row 4
subplot(4,2,7)
plot(p(:,1),p(:,2),'k')
ylabel('force')
xlabel('time')
axis([0 10 0 2])

subplot(4,2,8)
plot(t(:,1),t(:,2),'k')
xlabel('time')
axis([0 10 0 2])
suptitle('Drag force for flow over an in-line oscillating cylinder')
hold off
