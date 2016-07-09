%testing interpolation for pressure (one node is a neuman BC)
function [f, Ainv] = interpolateP(x, y, q)
a11 = 1;  a12 = x(1);  a13 = y(1); a14 = x(1)*y(1);
a21 = 1;  a22 = x(2);  a23 = y(2); a24 = x(2)*y(2);
a31 = 1;  a32 = x(3);  a33 = y(3); a34 = x(3)*y(3);
a41 = 1;  a42 = x(4);  a43 = y(4); a44 = x(4)*y(4);
% if abs(q(1))>35 
%     a11=0;
%     a12 = 1;
%     a13 = 1;
%     a14 = x(1)+y(1);
% end
% if abs(q(2))>35
    a21=0;
    a22 = 1;
    a23 = 1;
    a24 = x(2)+y(2);
% end
% if abs(q(3))>35
%     a31=0;
%     a32 = 1;
%     a33 = 1;
%     a34 = x(3)+y(3);
% end
% if abs(q(4))>35
    a41=0;
    a42 = 1;
    a43 = 1;
    a44 = x(4)+y(4);
% end
A = [a11 a12 a13 a14;
     a21 a22 a23 a24;
     a31 a32 a33 a34;
     a41 a42 a43 a44];

Ainv = inv(A);

a = Ainv*q';

f= @(X,Y) a(1) + a(2)*X + a(3)*Y + a(4)*X*Y;