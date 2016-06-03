function [f, Ainv] = interpolateP(x, y, q)
a11 = 1;  a12 = x(1);  a13 = y(1); a14 = x(1)*y(1);
a21 = 1;  a22 = x(2);  a23 = y(2); a24 = x(2)*y(2);
a31 = 1;  a32 = x(3);  a33 = y(3); a34 = x(3)*y(3);
a41 = 1;  a42 = x(4);  a43 = y(4); a44 = x(4)*y(4);
if q(1)==0
    a11=0;
    a12 = 1;
    a13 = 1;
    a14 = x(1)+y(1);
end
if q(2)==0
    a21=0;
    a22 = 1;
    a23 = 1;
    a24 = x(2)+y(2);
end
if q(3)==0
    a31=0;
    a32 = 1;
    a33 = 1;
    a34 = x(3)+y(3);
end
if q(4)==0
    a41=0;
    a42 = 1;
    a43 = 1;
    a44 = x(4)+y(4);
end
A = [a11 a12 a13 a14;
     a21 a22 a23 a24;
     a31 a32 a33 a34;
     a41 a42 a43 a44];

Ainv = inv(A);

a = Ainv*q';

f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;