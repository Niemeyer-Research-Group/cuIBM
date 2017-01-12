%look at divergence
clc
clear
close all
path = strcat('/scratch/src/cuIBM/validation/luo/test/Divergence.csv');
delim = '\t';
M = dlmread(path,delim,1,0);
%-v_north  dy_north  -u_east  dx_east  v_south  dy_south  u_west  dx_west  x  y
north = M(:,1)./M(:,2);
east = M(:,3)./M(:,4);
south = M(:,5)./M(:,6);
west = M(:,7)./M(:,8);
dy = M(:,2);
dx = M(:,4);
x = M(:,9);
y = M(:,10);
sum = north + east + south + west;
index = 1:length(north);
% for index = 1:length(north)
%     close
plot3(x(index),y(index),sum(index),'ok'), hold on
% plot3(x(index),y(index)+dy(index)/2,north(index),'xr')
% plot3(x(index)+dx(index)/2,y(index),east(index),'sb')
% plot3(x(index),y(index)-dy(index)/2,south(index),'^r')
% plot3(x(index)-dx(index)/2,y(index),west(index),'b+')

% plot cylinder
theta = linspace(0,pi*2,length(M));
Z = zeros(2,length(M));
X = Z;
Y = Z;
X(1,:) = cos(theta)*0.5;
X(2,:) = cos(theta)*0.5;
Y(1,:) = sin(theta)*0.5;
Y(2,:) = sin(theta)*0.5;
Z(1,:) = round(max(sum)+1);
Z(2,:) = round(min(sum)-1);
surf(X,Y,Z,'FaceColor','blue','EdgeColor','none'), hold on %plot cylinder surface
surf([-1 -1; 1 1], [-1 1; -1 1], [0 0; 0 0], 'FaceColor','blue','EdgeColor','black')
colormap hsv
alpha(0.4)



axis([min(x(index))-0.03, max(x(index))+0.03, min(y(index))-0.03, max(y(index))+0.03, min([sum(index); north(index); south(index); west(index); east(index)]), max([sum(index); north(index); south(index); west(index); east(index)])])
legend('sum','North','East','South','West')
xlabel('x')
ylabel('y')
% end