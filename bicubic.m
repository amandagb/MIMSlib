%=====================================================================
function f = bicubic(x)
%% Script: f = bicubic(x)
% Description: Function is copied from matlab's imresize function an
% outputs a bicubic function of x. The function is piecewise constant 
% INPUTS ----------------------------------------------------------------
% 
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author              E-mail                      Version
%  17 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%     Reads in data generated on 20141016 and creates a probability image
%     for the entire image (uses 5 classes: bkgnd, bld, EB, hemo, tsu
%  19 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1

% See Keys, "Cubic Convolution Interpolation for Digital Image
% Processing," IEEE Transactions on Acoustics, Speech, and Signal
% Processing, Vol. ASSP-29, No. 6, December 1981, p. 1155.

absx = abs(x);
absx2 = absx.^2;
absx3 = absx.^3;

f = (1.5*absx3 - 2.5*absx2 + 1) .* (absx <= 1) + ...
                (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) .* ...
                ((1 < absx) & (absx <= 2));
%---------------------------------------------------------------------