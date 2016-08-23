%% hybrid node velocity interpolation X 3D
clc
clear
close all
M = dlmread('/scratch/src/cuIBM/validation/osc/gh/interp_testX.csv','\t',1,0); %start on second row to avoid headers
% 1     2       3       4       5       6       7       8       9       10      11  12  13  14  15  16  17  18  19  20  21  22  23  24
% BN_X1	BN_Y1   BN_X2	BN_Y2	GN_X    GN_Y	BI_X	BI_Y    IP_X	IP_Y	x1	x2	x3  x4	y1	y2	y3	y4	q1	q2	q3	q4	GN_Uip_u
X = zeros(1,7);
Y = zeros(1,7);
Z = zeros(1,7);
for i =38%1:length(M)
%     close
    X(1) = M(i,5); %ghost node
    X(2) = M(i,7); %body intercept
    X(3) = M(i,11); %corner1
    X(4) = M(i,12); %corner2
    X(5) = M(i,13); %corner3
    X(6) = M(i,14); %corner4
    x1 = [M(i,1),M(i,3)]; %body
    Y(1) = M(i,6);
    Y(2) = M(i,8);
    Y(3) = M(i,15);
    Y(4) = M(i,16);
    Y(5) = M(i,17);
    Y(6) = M(i,18);
    y1 = [M(i,2),M(i,4)];
    Z(1) = M(i,23);
    Z(2) = -1.0; %should be set to body velocity
    Z(3) = M(i,19);
    Z(4) = M(i,20);
    Z(5) = M(i,21);
    Z(6) = M(i,22);
    z1 = [0,0];

%     Q = interpolate([X(3) X(4) X(5) X(6)], [Y(3) Y(4) Y(5) Y(6)], [Z(3) Z(4) Z(5) Z(6)]);
%     ip_u = Q(X(7),Y(7));
    
    scatter3(X(1),Y(1),Z(1),'ks'), hold on %ghost node
    scatter3(X(2),Y(2),Z(2),'ko') %body intercept
    plot3([X(3:4) X(6) X(5) X(3)], [Y(3:4) Y(6) Y(5) Y(3)], [Z(3:4) Z(6) Z(5) Z(3)],'kx-'); %interp corners
end
    
axis square
legend('Ghost node', 'Body Intercept', 'Corners')
% end
%% hybrid node velocity interpolation X 2d
clc
clear
close all
M = dlmread('/scratch/src/cuIBM/validation/osc/gh/interp_testX.csv','\t',1,0); %start on second row to avoid headers
% M = dlmread('/scratch/src/cuIBM/validation/cylinder/Re40/interp_testY.csv','\t',1,0); %start on second row to avoid headers
% 1     2       3       4       5       6       7       8       9       10      11  12  13  14  15  16  17  18  19  20  21  22  23  24
% BN_X1	BN_Y1   BN_X2	BN_Y2	GN_X    GN_Y	BI_X	BI_Y    IP_X	IP_Y	x1	x2	x3  x4	y1	y2	y3	y4	q1	q2	q3	q4	GN_Uip_u
X = zeros(1,7);
Y = zeros(1,7);
Z = zeros(1,7);
hold on
x = linspace(0,2*pi, 360);
plot (0.5*cos(x),0.5*sin(x),'g');
for i =1:length(M)
    X(1) = M(i,5); %ghost node
    X(2) = M(i,7); %body intercept
    X(3) = M(i,11); %corner1
    X(4) = M(i,12); %corner2
    X(5) = M(i,13); %corner3
    X(6) = M(i,14); %corner4
    X(7) = M(i,9); %image piont
    Y(1) = M(i,6);
    Y(2) = M(i,8);
    Y(3) = M(i,15);
    Y(4) = M(i,16);
    Y(5) = M(i,17);
    Y(6) = M(i,18);
    Y(7) = M(i,10);
    Z(1) = M(i,23);
    Z(2) = 0;
    Z(3) = M(i,19);
    Z(4) = M(i,20);
    Z(5) = M(i,21);
    Z(6) = M(i,22);
    Z(7) = M(i,24);
    z1 = [0,0];
    plot(X(1),Y(1),'ks')                            %Ghost node
    plot(X(2),Y(2),'ko')                            %body intercept
    plot(X(7),Y(7),'bx')                            %image point
    plot(X(3:6),Y(3:6),'rd')                        %interpolation corners
    plot(M(i,1),M(i,2),'rs',M(i,3),M(i,4),'rs')     %body nodes
    plot([X(1) X(7)], [Y(1) Y(7)], 'k-')            % line between ghost node and image point
    linecolor = rand(1,3);
    for j=3:6
        plot([X(j) X(7)], [Y(j) Y(7)], 'color', linecolor) % line between corner and image point
    end
end
hold off 
figure(1)
set(gca,'Color',[0.8 0.8 0.8]);
legend('Body','Ghost Node','Boundary Intercept','Image Point','Interpolation corner')
axis square