function [c] = sastatescalefactors(varargin)
%% Script: [c] = sastatescalefactors(varargin)
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  04 Mar 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

%=====================================================
% Simulated Annealing Parameters
%=====================================================
tx_scalefactor = 1;
ty_scalefactor = 1;
rot_scalefactor = 180/pi;
sx_scalefactor = 10;
sy_scalefactor = 10;
k_scalefactor = 50;
c = [tx_scalefactor,ty_scalefactor,rot_scalefactor,...
    sx_scalefactor,sy_scalefactor,k_scalefactor];
end