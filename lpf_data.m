function filtdata = lpf_data(data,filter_order,cutoff,varargin)
%% Script:  
% Description:  
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% data:  
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  24 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 1e3;
end

if ~exist('filter_order'); filter_order = 10; end
if ~exist('cutoff'); cutoff = 0.07; end

data = data(:);
data = [data(1:filter_order/2-1);data;data(end-filter_order/2+1:end)];

% Low pass filtering
h = fir1(filter_order,cutoff);
filtdata = filter(h,1,data);
filtdata = filtdata(filter_order:end);
end


