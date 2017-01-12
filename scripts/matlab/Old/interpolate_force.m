clc
clear
close all

%% force pressure interpolation 2D
clc
clear
close all
M = dlmread('/scratch/src/cuIBM/validation/luo/test/interp_test_force.csv','\t',1,0); %start on second row to avoid headers
% 1     2   3   4   5   6   7   8   9   10  11  12  13 14    15    16    17   18   19  20  21   22   23
% x1    x2  x3  x4  y1  y2  y3  y4	q1	q2	q3	q4	p  bnx1  bny1  bnx2  bny2 px   py  px2 py2  px3  py3
for i =1:length(M)
    linecolor = rand(1,3);
    plot(M(i,18),M(i,19),'ro'), hold on %body intercept
    plot(M(i,20),M(i,21),'kx') %cloeset node to bi
    plot(M(i,22),M(i,23),'gs') %image point
    for j=1:4
        plot(M(i,j),M(i,j+4),'rd')
%         plot([M(i,j) M(i,18)], [M(i,j+4) M(i,19)], 'color', linecolor)
        plot([M(i,j) M(i,22)], [M(i,j+4) M(i,23)], 'color', linecolor)
    end
    plot([M(i,22) M(i,18)], [M(i,23) M(i,19)],'-b')% ip to bi
%     plot([M(i,20) M(i,18)], [M(i,21) M(i,19)])% close to bi
    
end
plot(M(:,14),M(:,15),'-ko')
set(gca,'Color',[0.8 0.8 0.8]);
legend('Body Intercept', 'Closest node X', 'image point', 'Interpolation Corner')
xlabel('x')
ylabel('y')
axis square

%% force pressure interpolation 3D
clc
clear
close all
M = dlmread('/scratch/src/cuIBM/validation/luo/test/interp_test_force.csv','\t',1,0); %start on second row to avoid headers
% 1     2   3   4   5   6   7   8   9   10  11  12  13 14    15    16    17   18   19  20  21   22   23
% x1    x2  x3  x4  y1  y2  y3  y4	q1	q2	q3	q4	p  bnx1  bny1  bnx2  bny2 px   py  px2 py2  px3  py3
for i =1:length(M)
    linecolor = rand(1,3);
    scatter3(M(i,18),M(i,19),M(i,13),'ro'), hold on 
end
plot(M(:,14),M(:,15),'-ko')
set(gca,'Color',[0.8 0.8 0.8]);
legend('Body Intercept', 'Surface')
xlabel('x')
ylabel('y')
axis square