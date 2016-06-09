clear
clc
%close all
%test points if you don't have body point data
% x1 = 2; %lagrangian body point 1
% y1 = 10;
% x2 = 5; %lagrangian body point 2
% y2 = -1;
% x3 = 1;%ghost node
% y3 = 1;

%read data from cuIBM output, format is as follows
%x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
M = dlmread('/scratch/src/cuIBM/validation/luo/test/body_nodesX.csv','\t',1,0); %start on second row to avoid headers %change Y and X for different tests
%% compare the cuibm output data to an easier to manipulate version in matlab
for i =1:length(M)
index = i;
x1 = M(index,1);
y1 = M(index,2);
x2 = M(index,3);
y2 = M(index,4);
x3 = M(index,5);
y3 = M(index,6);

p = sqrt((x1-x2)^2 + (y1-y2)^2); %distance from 1 to 2
o = sqrt((x1-x3)^2 + (y1-y3)^2); %distance from 3 to 1
b = sqrt((x3-x2)^2 + (y3-y2)^2); %distance from 2 to 3

theta_321 = acos((p^2+b^2-o^2)/(2*p*b)); %angle format: start point,vertex,endpoint

a = sin(theta_321)*b;%distance from 3 to 4 
%instead of finding a then doing a/tan(theta) I think you can skip one of
%these steps

y4 = y2 + (y1-y2)/p * a/tan(theta_321); %body intercept point, forms a right angle with lines 1-2 and 4-3
x4 = x2 + (x1-x2)/p * a/tan(theta_321);
y5 = y4+(y4-y3); %image point, mirror the ghost node accross line 1-2
x5 = x4+(x4-x3);

xp1 = [x1,x2]; %make lines to plot in matlab
yp1 = [y1,y2];
xp2 = [x2,x3];
yp2 = [y2,y3];
xp3 = [x3,x1];
yp3 = [y3,y1];
xp4 = [x3,x5];
yp4 = [y3,y5];
plot(x1,y1,'rs',x3,y3,'ks',x4,y4,'ko',x5,y5,'k^',M(index,7),M(index,8),'gx',M(index,9),M(index,10),'rx',x2,y2,'rs'), hold on
plot(xp1,yp1,'k',xp2,yp2,'k',xp3,yp3,'k',xp4,yp4,'k'),hold on

%check that we made a right angle, print if it isn't a right angle
p = sqrt((x2-x4)^2 + (y2-y4)^2); %distance from 2 to 4
o = sqrt((x2-x3)^2 + (y2-y3)^2); %distance from 2 to 3
b = sqrt((x3-x4)^2 + (y3-y4)^2); %distance from 3 to 4

theta_243 = acos((p^2+b^2-o^2)/(2*p*b));

if abs(rad2deg(theta_243)-90) > 0.0001
    fprintf('point i%d does not make a right angle: %f\n', i, rad2deg(theta_243))
end
if (x5 - M(i,10) > 0.00001 && y5 - M(i,10) > 0.00001)
    fprintf('cuIBM and matlab dont match at index i%d   (%03f,%03f)\n',i,x5,y5)
    plot(x5,y5,'^')
end
end
legend('Lagrangian Body Tracker','Ghost Node','Body Intercept','Image Point','C++ BI','C++ IP')
axis square

%% plot cuibm output
figure(2)
plot(M(:,5),M(:,6),'rs',M(:,7),M(:,8),'ko',M(:,9),M(:,10),'kx'), hold on
plot(0.5*sin(0:0.01:2*pi),0.5*cos(0:0.01:2*pi),'k')
legend('Ghost Node','Boundary Intercept','Image Point')
%set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','g') %%set fill color
axis square
axis([-0.6,0.6,-0.6,0.6])





