function [A,B,iB] = alpha2tmat(alpha,M,N)
% [ysize,xsize] = size(I) => ysize = # of rows, xsize = # of columns

tx = alpha(1);
ty = alpha(2);
theta = alpha(3);
if length(alpha) > 3
  sx = alpha(4);
  sy = alpha(5);
else
  sx = 1;
  sy = 1;
end
if length(alpha) > 5
  k = alpha(6);
end

Tt = [1,0,tx; 0,1,ty; 0,0,1];
Tr = [cos(theta),-sin(theta),0; sin(theta),cos(theta),0; 0,0,1];
Tg = [1,k,0;0,1,0;0,0,1];
Ts = [sx,0,0; 0,sy,0; 0,0,1];
A = Tt*Tr*Tg*Ts;

x_center = N/2;
y_center = M/2;
B = A';
B(3) = B(3) - B(1)*x_center - B(2)*y_center + x_center;
B(6) = B(6) - B(4)*x_center - B(5)*y_center + y_center;

iB = inv(B);
iB(7:9) = [0;0;1];
end