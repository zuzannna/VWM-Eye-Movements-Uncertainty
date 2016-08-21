function [x_c,y_c] = zuzaspenthingtoday(angle)

%constant = 2*x_f - 2*tan(theta)*y_f-d^2+x_f^2+y_f^2;

%x_c = -sqrt(-4*constant)./2;

%y_c = tan(theta)*x_c;

targetRadius=10;
%angle=50;

x_c = targetRadius * cos(pi*angle/180);
y_c = targetRadius * sin(pi*angle/180);