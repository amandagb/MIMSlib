function [mu,phiMI] = fmin_MI(Fin,Min,varargin)
%% Script:
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% Fin:  fixed image
% Min:  moving image
% varargin - 'PropertyName','PropertyValue'
%   'mu0':    1x6 vector with initial transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  24 Oct 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

global M nbins klog_cell pF
M = Min;
F = Fin;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('mu0',PropertyNames)
  mu0 = PropertyVal{strmatch('mu0',PropertyNames)};
else
  mu0 = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk
end

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 50;
end

% Initialize fixed image
[Fm,Fn] = size(F);
LF = nbins; % Total number of histogram used to approximate pdf of F
Frng = max(F(:)) - min(F(:));
bin_rngF = Frng/(LF-1);
Findex = (F-min(F(:)))./bin_rngF;
% Compute 0-spline histogram of fixed image
pF = zeros(1,LF);
klog_cell = cell(LF,1);
[klog_cell{1:end}] = deal(zeros(Fm,Fn));
for k = 1:LF
  b0arg = (k-1) - Findex;
  klog_cell{k} = (b0arg >= -1/2 & b0arg < 1/2); %find(b0arg <= 0 & b0arg > -1);
  pF(k) = sum(klog_cell{k}(:));
end
normfac = (sum(pF));
pF = pF./normfac;
optfmin = optimset('PlotFcns',{@optimplotx,@optimplotfunccount,@optimplotfval,...
  @optimplotstepsize,@optimplotfirstorderopt});
[mu,phiMI] = fminunc(@computeMI,mu0,optfmin);

  function SMI = computeMI(mu); % would need M, nbins, klog_cell, pF, LF to be global
    MI = twoimgMI([],M,'mu',mu,'klog_cell',klog_cell,'pA',pF,'nbins',nbins);
    SMI = -MI;
  end
end


































