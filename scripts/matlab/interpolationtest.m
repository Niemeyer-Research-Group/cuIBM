clc
clear
close all

x1 = 10;
y1 = 10;
q1 = 0;

x2 = 11;
y2 = 10;
q2 = 1;

x3 = 10;
y3 = 11;
q3 = 2;

x4 = 10.5;
y4 = 10.5;
q4 = 1;

a0 = q1;
a1 = 10;
a2 = 10;
a3 = 10;
a30 = 30;
a20 = 30;
a10 = 30;
flag = true;
tol = 0.001;
count = 0;
alpha = 0.75;
while flag
    a0 = ((q1     -a1*x1  - a2*y1 - a3*x1*y1)        )* alpha + (1-alpha)*a0;
    a1 = ((q2 -a0         - a2*y2 - a3*x2*y2)/x2     )* alpha + (1-alpha)*a1;
    a2 = ((q3 -a0 - a1*x3         - a3*x3*y3)/y3     )* alpha + (1-alpha)*a2;
    a3 = ((q4 -a0 - a1*x4 - a2*y4           )/(y4*x4))* alpha + (1-alpha)*a3;

    if (abs(a3 - a30) < tol && abs(a2 - a20) < tol && (a1 - a10) < tol || count > 10000)
        flag = false;
        count
    end
    count = count +1;
    a30 = a3;
    a20 = a2;
    a10 = a1;
end

Q= @(x,y) a0 + a1*x + a2*y + a3*x*y;
n = 10;
x = linspace(x1,x2,n);
y = linspace(y1,y3,n);
for i=1:length(x)
    for j = 1:length(y)
        z = Q(x(i),y(i)); 
        if (z >max([q1,q2,q3,q4]) || z < min([q1,q2,q3,q4]))
            fprintf('error q out of bounds\n');
        end
        scatter3(x(i),y(j),Q(x(i),y(j)),'ro'), hold on
    end
end

X = [x1 x2 x3 x4];
Y = [y1 y2 y3 y4];
Z = [q1 q2 q3 q4];
scatter3(X,Y,Z,'kx')
xlabel('x')
ylabel('y')
