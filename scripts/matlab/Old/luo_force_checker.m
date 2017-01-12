%This script is used to debug the luo force calculation kernels. It needs the files interp_test_force.csv and inter_test_force_dudn.csv 
%these are generated with the functions testForce_p() and testForce_dudn() currently commented out in the function luoIBM/calculateForce.inl luoForce()

close
clear
clc

%% pressure
figure
path = '/scratch/src/cuIBM/validation/cylinder/Re40/interp_test_force.csv';
delim = '\t';
M = dlmread(path,delim,1,0);


hold on
for i=1:length(M)
    scatter3(M(i,1),M(i,2),M(i,3))
end
hold off
xlabel('x')
ylabel('y')
zlabel('p')

%% velocity
figure
path = '/scratch/src/cuIBM/validation/cylinder/Re40/interp_test_force_dudn.csv';
delim = '\t';
M = dlmread(path,delim,1,0);

hold on
for i=1:length(M)
    scatter3(M(i,1),M(i,2),M(i,5),'r')
%     scatter3(M(i,1),M(i,2),M(i,3),'r')
end
hold off
xlabel('x')
ylabel('y')
% zlabel('dudn')
zlabel('pressure * area * n1')
