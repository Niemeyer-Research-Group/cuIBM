%the inside pressure wasn't working very well so this checks it
clc
clear
close all

%import data
body = dlmread('/scratch/src/cuIBM/validation/luo/test/body.csv','\t');
% x, y
inside = dlmread('/scratch/src/cuIBM/validation/luo/test/insidePressure.csv','\t');
% x, y, p
outside = dlmread('/scratch/src/cuIBM/validation/luo/test/outsidePressure.csv','\t');
%make cylinder
% x, y, p
theta = linspace(0,pi*2,length(body));
z = zeros(2,length(body));
x = z;
y = z;
x(1,:) = cos(theta)*0.5;
x(2,:) = cos(theta)*0.5;
y(1,:) = sin(theta)*0.5;
y(2,:) = sin(theta)*0.5;
z(1,:) = round(max(outside(:,3))+1);
z(2,:) = round(min(outside(:,3))-1);
surf(x,y,z,'FaceColor','blue','EdgeColor','none'), hold on %plot cylinder surface
colormap hsv
alpha(0.4)
plot3(outside(:,1),outside(:,2),outside(:,3),'ko')%plot outside nodes
plot3(inside(:,1),inside(:,2),inside(:,3),'ro')%plot inside nodes
plot(body(:,1),body(:,2))%plot cylinder
legend('Surface','Outside', 'Inside')
xlabel('x')
ylabel('y')
zlabel('Pressure')
axis([-1,1,-1,1])

%% prove that values that are too far away are the incorrect ones
% close
% count = 1;
% for i=1:length(inside)
%     distance = sqrt(inside(i,1)^2+inside(i,2)^2);
%     if 0.5-distance >0.01
%         INSIDE(count,1) = inside(i,1);
%         INSIDE(count,2) = inside(i,2);
%         INSIDE(count,3) = inside(i,3);
%         count = count + 1;
%     end
% end
% 
% surf(x,y,z,'FaceColor','blue','EdgeColor','none'), hold on %plot cylinder surface
% colormap hsv
% alpha(0.4)
% plot3(outside(:,1),outside(:,2),outside(:,3),'ko')%plot outside nodes
% plot3(INSIDE(:,1),INSIDE(:,2),INSIDE(:,3),'ro')%plot inside nodes
% plot(body(:,1),body(:,2))%plot cylinder
% legend('Surface','Outside', 'Inside')
% xlabel('x')
% ylabel('y')
% zlabel('Pressure')
% axis([-1,1,-1,1])
% 
% %% plot full pressure body
% clear all
% close all
% clc
% inside = dlmread('/scratch/src/cuIBM/validation/luo/test/3dPressure.csv','\t');
% plot3(inside(:,1),inside(:,2),inside(:,3),'ko')
% xlabel('x')
% ylabel('y')
% zlabel('Pressure')
% axis([-1,1,-1,1])