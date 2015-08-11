function [thresh_ind,min_changes,interm_thresh,low_ind,varargout] = find_change_thresh(dvec,varargin);

%% Script: [thresh_ind,min_changes,interm_thresh,low_ind,varargout] = find_change_thresh(dvec,varargin);
% Description: Given a vector of data (dvec) with notably different
% foreground and background, the script finds the appropriate threshold
% which divides the data. The function can find the threshold starting from
% the top of the data vector values or from the bottom. Indicating
% 'starting_point' as 't' (meaning top) will give a interm_thresh greather
% than when 'starting_point' is indicated as 'b' (meaning start searching
% from the bottom values of the data vector).
%      ______________
%      |             |
%      |             |
% _____|             |_____
%   
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'starting_point': string which inidiates whether to begin the search
%   from the top (maximum) to the bottom (minimum) of the vector => 't'; or
%   from the bottom to the top of the sorted vector values => 'b';
%   'nchanges': number of times the vector is expected to change between
%   low and high values [DEFAULT = 2 => low to high to low]
% OUTPUTS ---------------------------------------------------------------
% la0info: structure which contains information about laser off periods
%   la0ind:
%   meanla0ind:
%   la0peak:
%   meanla0peak:
%
%  Date           Author            E-mail                      Version
%  24 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('start',PropertyNames)
  start_str = PropertyVal{strmatch('start',PropertyNames)};
  if strmatch('t',start_str)
    start_top = 1;
  else 
    start_top = 0;
  end
else
  start_top = 0;
end

if strmatch('nchanges',PropertyNames)
  reqchngs = PropertyVal{strmatch('nchanges',PropertyNames)};
else
  reqchngs = 2;
end

Ddvec = diff(dvec);
sortdvec = sort(dvec);
min_changes = 100;
thresh_ind = [0,0];
interm_thresh = 0;
binaryvec = zeros(size(dvec));

if start_top
  [v,threshi] = min(abs(sortdvec-mean(dvec)));
  thresh = sortdvec(threshi);
else
  sortdvec = flipud(sortdvec(:));
  threshi = length(sortdvec) - 10;
  thresh = sortdvec(threshi);
end

highlowfail = 1;
while highlowfail && threshi > 1
  binaryvec = dvec > thresh;
  chngs = find(diff(binaryvec) ~= 0);
  if length(chngs) == 1
    highlowfail = 0;
    low_ind = find(binaryvec == 0);
    if binaryvec(chngs(1)) == 0
      thresh_ind = [chngs,0];
    else
      thresh_ind = [0,chngs];
    end
  elseif length(chngs) == 2 && binaryvec(chngs(1)) == 0
    highlowfail = 0;
    low_ind = find(binaryvec == 0);
    thresh_ind = chngs;
  else
    threshi = threshi - 1;
    thresh = sortdvec(threshi);
  end
  
  if length(chngs) < min_changes && rem(length(chngs),2) == 0 ...
          || length(chngs) == 1
    min_changes = length(chngs);
    interm_thresh = sortdvec(threshi);
  end
end
if rem(length(chngs),2) == 0;
  min_changes = min(length(chngs),min_changes);
  interm_thresh = thresh;
end
end