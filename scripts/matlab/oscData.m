clear
close all
clc

%pixels for data points
x = [398 540 682 824 965 1108 1250 1392 1534 1675 1818];
y = [405 550 705 776 776 698 541 392 320 319 405];

%real data points
xx = zeros(length(x),1);
yy = xx;

%pixels for grid
X=[1818 398];
Y=[832 263];

%real grid
XX =[1 0];
YY = [-2 2];

slopex = (XX(1)-XX(2))/(X(1)-X(2));
slopey = (YY(1)-YY(2))/(Y(1)-Y(2));
for i=1:length(x)
    xx(i) = slopex*(x(i)-X(2));
    yy(i) = slopey*(y(i)-Y(1)) - 2;
end



