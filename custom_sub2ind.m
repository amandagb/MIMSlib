function linind = custom_sub2ind(SZX,x)
%% Script: linind = custom_sub2ind(SZX,x)
% Description: Extension of sub2ind function
% Example:
% INPUTS ----------------------------------------------------------------
% 
% OUTPUTS ---------------------------------------------------------------
% 
%  Date           Author            E-mail                      Version
%  7  Oct 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

if length(SZX) == 2
  linind = sub2ind(SZX,x(1),x(2));
elseif length(SZX) == 3
  linind = sub2ind(SZX,x(1),x(2),x(3));
elseif length(SZX) == 4
  linind = sub2ind(SZX,x(1),x(2),x(3),x(4));
elseif length(SZX) == 5
  linind = sub2ind(SZX,x(1),x(2),x(3),x(4),x(5));
elseif length(SZX) == 6
  linind = sub2ind(SZX,x(1),x(2),x(3),x(4),x(5),x(6));
end
