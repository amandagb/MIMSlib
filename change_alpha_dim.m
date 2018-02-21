function [alpha,A,alphainv,Ainv] = change_alpha_dim(alpha,oldsize,newsize)
%% Script: [alpha,A,alphainv,Ainv] = change_alpha_dim(alpha,oldsize,newsize)
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
%  21 Nov  2016   Amanda Balderrama amanda.gaudreau@gmail.com     1

sizeScale = newsize./oldsize;

[mu,A,~,~] = find_inv_alpha(alpha,oldsize(1),oldsize(2));
Tupsamp = [sizeScale(2),0,0; 0,sizeScale(1),0; 0,0,1];
Tdownsamp = [1/sizeScale(2),0,0; 0,1/sizeScale(1),0; 0,0,1];
[munew,~,~,~] = find_inv_alpha(Tupsamp*A*Tdownsamp,newsize(1),newsize(2));
% munew(4) = mu(4); %munew(4)/sizeScale(2); %
% munew(5) = mu(5); %munew(5)/sizeScale(1); %

[alpha,A,alphainv,Ainv] = find_inv_alpha(munew,newsize(1),newsize(2));
end