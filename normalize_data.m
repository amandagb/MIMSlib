function [dout,th] = normalize_data(din,thresh,varargin)
%% Script: [normd,th] = normalize_data(din,thresh,varargin)
% Description: Normalizes that data in the 2D matrix din to fall between 0
% and 1
% Example:  dNmat = normalize_data(dmat,[]);
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'norm': binary indicating whether to normalize (1) or not (0) the data
%   to the range 0 to 1 [DEFAULT = 1]
% OUTPUTS ---------------------------------------------------------------
% 
%
%  Date           Author            E-mail                      Version
%  22 June 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

if ~exist('thresh'); thresh = []; end;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

nchannels = size(din,3);
dout = zeros(size(din));

if strmatch('norm',PropertyNames)
  n = PropertyVal{strmatch('norm',PropertyNames)};
else 
  n = 1;
end

for c = 1:nchannels
  dch = din(:,:,c);
  if ~thresh
    th = [];
  elseif isempty(thresh)
    th = auto_thresh(dch,'auto',[]);
    dch(dch <= th(1)) = th(1);
    dch(dch >= th(2)) = th(2);
  else
    th = thresh;%auto_thresh(din,[thresh],[]);
    dch(dch <= th(1)) = th(1);
    dch(dch >= th(2)) = th(2);
  end
  
  if n
    normd = (dch - min(dch(:)))./(max(dch(:))-min(dch(:)));
  else
    normd = dch;
  end
  dout(:,:,c) = normd;
end
end