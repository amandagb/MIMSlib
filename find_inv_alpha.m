function [alpha,A,alphainv,Ainv] = find_inv_alpha(alpha,M,N)
%% Script: [alphainv] = find_inv_alpha(alpha)
% Description: [tx=0; ty=0; theta=0; sx=1; sy=1; k=0]
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% alpha:  1x6 vector with the transformation parameters OR
%         3x3 matrix transformation matrix
%         [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew]
%         NOTE: Positive displacement is to the left, up, and counterclockwise,
%           respectively.
% varargin - 'PropertyName','PropertyValue'
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  6  Dec 2012    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  17 Oct 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     2

y_center = (M+1)/2;
x_center = (N+1)/2;
if size(alpha) == [3,3]
  if sum(alpha(3,:) == [0,0,1]) == 1  % input matrix is of the form M
    Mmatlab = alpha;
    A = Mmatlab';
    A(1,3) = Mmatlab(3) + Mmatlab(1)*x_center + Mmatlab(2)*y_center - x_center;
    A(2,3) = Mmatlab(6) + Mmatlab(4)*x_center + Mmatlab(5)*y_center - y_center;
  else
    A = alpha;
  end
else % This case is executed if the elements of the transform vector are entered as a 6 element vector
  [A,Mmatlab,iMmatlab] = alpha2tmat(alpha,M,N);
end

a1 = A(1,1);
a2 = A(1,2);
a3 = A(1,3);
a4 = A(2,1);
a5 = A(2,2);
a6 = A(2,3);
tx = a3;   ty = a6;   rot = atan(a4/a1);
sx = a1/cos(rot);
sk = (a5*sin(rot) + a2*cos(rot))/(a5*cos(rot) - a2*sin(rot));
sy = a5/(sk*sin(rot) + cos(rot));
alpha = [tx,ty,rot,sx,sy,sk];

Ainv = inv(A);
b1 = Ainv(1,1);
b2 = Ainv(1,2);
b3 = Ainv(1,3);
b4 = Ainv(2,1);
b5 = Ainv(2,2);
b6 = Ainv(2,3);
txp = b3;   typ = b6;   rotp = atan(b4/b1);
sxp = b1/cos(rotp);
skp = (b5*sin(rotp) + b2*cos(rotp))/(b5*cos(rotp) - b2*sin(rotp));
syp = b5/(skp*sin(rotp) + cos(rotp));
alphainv = [txp,typ,rotp,sxp,syp,skp];
end