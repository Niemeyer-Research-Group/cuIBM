%% 3d pressure interpolation
clc
clear
% close all
figure
M = dlmread('/scratch/src/cuIBM/validation/osc/flow/interp_testP.csv','\t',1,0); %start on second row to avoid headers
% M = dlmread('/scratch/src/cuIBM/validation/cylinder/Re40/interp_testP.csv','\t',1,0); %start on second row to avoid headers
% 1     2       3       4       5       6       7       8       9       10      11  12  13  14  15  16  17  18  19  20  21  22  23 24 25 26 27 28
% BN_X1	BN_Y1   BN_X2	BN_Y2	GN_X    GN_Y	BI_X	BI_Y    IP_X	IP_Y	x1	x2	x3  x4	y1	y2	y3	y4	q1	q2	q3	q4	p* a0 a1 a2 a3 BI_p
X = zeros(1,7);
Y = zeros(1,7);
Z = zeros(1,7);
hold on
for i = 45%1:length(M)
    X(1) = M(i,5); %ghost node
    X(2) = M(i,7); %body intercept
    X(3) = M(i,11); %corner1
    X(4) = M(i,12); %corner2
    X(5) = M(i,13); %corner3
    X(6) = M(i,14); %corner4
    X(7) = M(i,9); %image piont
    x1 = [M(i,1),M(i,3)]; %body
    Y(1) = M(i,6);
    Y(2) = M(i,8);
    Y(3) = M(i,15);
    Y(4) = M(i,16);
    Y(5) = M(i,17);
    Y(6) = M(i,18);
    Y(7) = M(i,10);
    y1 = [M(i,2),M(i,4)];
    Z(1) = M(i,23);
    Z(2) = 0;
    Z(3) = M(i,19);
    Z(4) = M(i,20);
    Z(5) = M(i,21);
    Z(6) = M(i,22);
    z1 = [0,0];
%     matD = M(i,34)+M(i,29)+M(i,30)+M(i,31)+M(i,32)+M(i,33)
%     scatter3(M(i,7),M(i,8),-matD)
%     [Q, a] = interpolateP(X(3:6), Y(3:6), Z(3:6));
    Q= @(X,Y) M(i,24) + M(i,25)*X + M(i,26)*Y + M(i,27)*X*Y;
    
    aa = scatter3(X(1),Y(1),Z(1),'ks');%ghost node
    cc = scatter3(X(2),Y(2),M(i,28),'kx'); %body intercept
    xx = linspace(X(3),X(4), 10);
    yy = linspace(Y(3),Y(5),10);
    for j = 1:length(xx)
        for k = 1:length(yy)
            dd = scatter3(xx(j),yy(k),Q(xx(j),yy(k)));
        end
    end
    for j=3:6 %corners
        if abs(Z(j))<35
            ee = scatter3(X(j),Y(j),Z(j),'rd');
        end
    end
%     plot3([X(1) X(2)], [Y(1) Y(2)], [Z(1),Z(2)],'--') %line between gn and cpp ip
end
axis square
<<<<<<< HEAD
% legend('Ghost node','Image Point', 'Body Intercept', 'Corners')
legend([aa cc dd ee],'ghost','BI','field','corner')
=======
legend('Ghost node', 'Body Intercept')
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
xlabel('x')
ylabel('y')
zlabel('pressure')
% end

<<<<<<< HEAD
%% 2d pressure interpolation
clc
clear
% close all
figure
%M = dlmread('/scratch/src/cuIBM/validation/osc/flow/interp_testP.csv','\t',1,0); %start on second row to avoid headers
M = dlmread('/scratch/src/cuIBM/validation/cylinder/Re40/interp_testP.csv','\t',1,0); %start on second row to avoid headers

% 1     2       3       4       5       6       7       8       9       10      11  12  13  14  15  16  17  18  19  20  21  22  23
% BN_X1	BN_Y1   BN_X2	BN_Y2	GN_X    GN_Y	BI_X	BI_Y    IP_X	IP_Y	x1	x2	x3  x4	y1	y2	y3	y4	q1	q2	q3	q4	p*
X = zeros(1,7);
Y = zeros(1,7);
Z = zeros(1,7);
for i =1:length(M)
    X(1) = M(i,5); %ghost node
    X(2) = M(i,7); %body intercept
    X(3) = M(i,11); %corner1
    X(4) = M(i,12); %corner2
    X(5) = M(i,13); %corner3
    X(6) = M(i,14); %corner4
    X(7) = M(i,9); %image point
    x1 = [M(i,1),M(i,3)]; %body
    Y(1) = M(i,6);
    Y(2) = M(i,8);
    Y(3) = M(i,15);
    Y(4) = M(i,16);
    Y(5) = M(i,17);
    Y(6) = M(i,18);
    Y(7) = M(i,10);
    y1 = [M(i,2),M(i,4)];
    Z(1) = M(i,23);
    Z(2) = 0;
    Z(3) = M(i,19);
    Z(4) = M(i,20);
    Z(5) = M(i,21);
    Z(6) = M(i,22);
    z1 = [0,0];
    plot(X(1),Y(1),'ks'),hold on %Ghost node
    plot(X(2),Y(2),'ko'), hold on %body intercept
    plot(x1,y1,'g-'), hold on %line that is the body
    plot(X(7),Y(7),'bx'), hold on %image point
    plot(X(3:6),Y(3:6),'rd'), hold on %interpolation corners
    plot(M(i,1),M(i,2),'rs',M(i,3),M(i,4),'rs') %body nodes
    %plot([X(1) X(7)], [Y(1) Y(7)], 'k-') % line between ghost node and image point
    trouble_nodes = [45];
    if(any(i==trouble_nodes))
        linecolor = [1 0 0]; %if trouble spot set red else set blue
    else
        %linecolor = rand(1,3);
        linecolor = [0 0 1];
    end
    for j=3:6
       % plot([X(j) X(1)], [Y(j) Y(1)], 'color', linecolor) % line between corner and image point outside
       plot([X(j) X(7)], [Y(j) Y(7)], 'color', linecolor) % line between corner and image point inside
    end
end
    
set(gca,'Color',[0.8 0.8 0.8]);
legend('Ghost Node','Boundary Intercept','Body','Image Point','Interpolation corner')
axis square
=======
% %% 2d pressure interpolation
% clc
% clear
% close all
% M = dlmread('/scratch/src/cuIBM/validation/luo/test/interp_testP.csv','\t',1,0); %start on second row to avoid headers
% % 1     2       3       4       5       6       7       8       9       10      11  12  13  14  15  16  17  18  19  20  21  22  23
% % BN_X1	BN_Y1   BN_X2	BN_Y2	GN_X    GN_Y	BI_X	BI_Y    IP_X	IP_Y	x1	x2	x3  x4	y1	y2	y3	y4	q1	q2	q3	q4	p*
% X = zeros(1,7);
% Y = zeros(1,7);
% Z = zeros(1,7);
% for i =1:length(M)
%     X(1) = M(i,5); %ghost node
%     X(2) = M(i,7); %body intercept
%     X(3) = M(i,11); %corner1
%     X(4) = M(i,12); %corner2
%     X(5) = M(i,13); %corner3
%     X(6) = M(i,14); %corner4
%     X(7) = M(i,9); %image piont
%     x1 = [M(i,1),M(i,3)]; %body
%     Y(1) = M(i,6);
%     Y(2) = M(i,8);
%     Y(3) = M(i,15);
%     Y(4) = M(i,16);
%     Y(5) = M(i,17);
%     Y(6) = M(i,18);
%     Y(7) = M(i,10);
%     y1 = [M(i,2),M(i,4)];
%     Z(1) = M(i,23);
%     Z(2) = 0;
%     Z(3) = M(i,19);
%     Z(4) = M(i,20);
%     Z(5) = M(i,21);
%     Z(6) = M(i,22);
%     z1 = [0,0];
%     plot(X(1),Y(1),'ks'),hold on %Ghost node
%     plot(X(2),Y(2),'ko'), hold on %body intercept
%     plot(x1,y1,'g-'), hold on %line that is the body
%     plot(X(7),Y(7),'bx'), hold on %image point
%     plot(X(3:6),Y(3:6),'rd'), hold on %interpolation corners
%     plot(M(i,1),M(i,2),'rs',M(i,3),M(i,4),'rs') %body nodes
%     %plot([X(1) X(7)], [Y(1) Y(7)], 'k-') % line between ghost node and image point
%     linecolor = rand(1,3);
%     for j=3:6
%        % plot([X(j) X(1)], [Y(j) Y(1)], 'color', linecolor) % line between corner and image point outside
%        plot([X(j) X(7)], [Y(j) Y(7)], 'color', linecolor) % line between corner and image point inside
%     end
% end
%     
% figure(1)
% set(gca,'Color',[0.8 0.8 0.8]);
% legend('Ghost Node','Boundary Intercept','Body','Image Point','Interpolation corner')
% axis square
>>>>>>> parent of 1831b5e... luo method works for all reynolds numbers for the stationary cylinder
