function [thresh,pnt_ind,dth] = auto_thresh(d,x,pcORvalstr,varargin)
%% Script: [thresh,pnt_ind] = auto_thresh(d,x,nbins)
% Description: Backward search to find the intensity data that makes up the highest x-% of values
% Example:
% INPUTS ----------------------------------------------------------------
% d:  l x s matrix containing intensity values for the given element (l = #
% lines, s = # samples)
% x:  1 x 2 matrix where [%low_intensities, %high_intensities] are to be
% eliminated to define thresholds min and max thresholds. This can also be
% indicated as a string 'auto' where the default is used
%   [DEFAULT = bottom 0.02% and top 0.5%]
% pcORvalstr:  string indicating whether x is a percentage ('perc') or
% values indicaing a threshold ('val') [DEFAULT = 'perc']
% varargin - 'PropertyName','PropertyValue'
%   • 'histmethod':  logical indicating whether the histogram approximation
%         method (cerca 2012) should be used or not [DEFAULT = 0 --- new
%         mthod]
%             IF 'histmethod',1 THEN pcorvalstr instead should represent
%             the number of bins to use in the histogram method
%
% OUTPUTS ---------------------------------------------------------------
% thresh: 1 x 2 matrix where [min_thresh, max_thresh] are indicated for
% which low_int-% of data is < min_thresh and high_int-% data is >
% max_thresh
% pnt_ind:  1 x 2 cell where each element contains a N x 2 matrix.
% first element indicates line (row) index of points >
% thresh
% dth: l x s matrix with values in the original matrix thresholded
%
%  Date           Author            E-mail                      Version
%  21 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  20 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  14 Nov  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%  23 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.01
%   Changed condition in line 99
%   3 Apr  2015   AG Balderrama     amanda.gaudreau@gmail.com     4
%   Changed method for estimating max and min threshold

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if ~exist('pcORvalstr') || isempty(pcORvalstr)
  pcORvalstr = 'perc';
elseif isnumeric(pcORvalstr)
  pcORvalstr = 'perc';
  binmethod = 1;
end

if strmatch('histmethod',PropertyNames)
  binmethod = PropertyVal{strmatch('histmethod',PropertyNames)};
else
  binmethod = 0; %CHANGED TO PNG from 'tif' 23 Feb 2015
end

if ~exist('x') || isempty(x)
  x = 'auto';
end

if strmatch('auto',x)
  x = [2e-2,5e-1];
end

d_vec = d(:);
if binmethod == 0
  thresh = quantile(d_vec,[x(1)*1e-2,1-x(2)*1e-2]);
  min_t = thresh(1); max_t = thresh(2);
  [l,s] = find(d <= min_t);
  pnt_ind{1} = [l,s];
  [l,s] = find(d >= max_t);
  pnt_ind{2} = [l,s];
else
  error = 0;  % Can be set to a non-zero value
  npnts = size(d,1)*size(d,2);
  N = x.*npnts/100;   % Translates percentage of data in number of points, N
  
  [bin_cnts,x_loc] = hist(d_vec,nbins);
  x_spacing = x_loc(2) - x_loc(1);
  half_space = x_spacing/2;
  
  %% Search for N(1) points in d < x(1)
  if N(1)
    cumsum_bins = cumsum(bin_cnts);
    ind_rng = find(cumsum_bins <= N(1));
    if isempty(ind_rng); ind_rng = 0; end
    closest_ind = ind_rng+1;
    min_t = 0;
    fine_bcnts = N(1)+1;
    last_fbcnts = 0;
    while ~min_t
      [fine_bcnts,fine_xloc] = hist(d_vec(d_vec < (x_loc(closest_ind(end))+half_space)),cumsum_bins(closest_ind(end)));
      cumsum_bins = cumsum(fine_bcnts);
      [v,i] = min(abs(N(1) - cumsum_bins));
      closest_xloc = fine_xloc(i);
      if v == N(1)
        min_t = 0.1;
      elseif length(fine_bcnts) > 1
        data_inbin = d_vec(d_vec >= min(d_vec(:)) ... closest_xloc - 0.5*(fine_xloc(2)-fine_xloc(1)) ...
          & d_vec <= closest_xloc + 0.5*(fine_xloc(2)-fine_xloc(1)));
        min_t = max(data_inbin);
        if ~min_t;
          min_t = 0.1;
        end
      elseif d_vec(d_vec < (x_loc(closest_ind(end))+half_space)) == 0
        min_t = x_loc(closest_ind(end));
      else
        min_t = fine_xloc;
      end
    end
    if min_t == 0.1
      min_t = min(d_vec);
    else
      [l,s] = find(d <= min_t);
      pnt_ind{1} = [l,s];
    end
  else
    min_t = min(d_vec);
  end
  
  %% Search for N(2) points in d > x(2)
  bkwrd_bcnts = fliplr(bin_cnts);
  bkwrd_xloc = fliplr(x_loc);
  cumsum_bins = cumsum(bkwrd_bcnts);
  ind_rng = find(cumsum_bins <= N(2));
  if isempty(ind_rng); ind_rng = 0; end
  closest_ind = ind_rng(end)+1;
  half_space = (bkwrd_xloc(1) - bkwrd_xloc(2))/2;
  bin_range = [bkwrd_xloc(closest_ind) - half_space,...
    max(d_vec)];
  dbn = d_vec(d_vec >= bin_range(1) & d_vec <= bin_range(2));
  nbns = length(dbn);
  pnt_cnt = 0;
  max_t = 0;
  while ~max_t
    [fine_bcnts,fine_xloc] = hist(dbn,nbns);
    closest_ind = find(cumsum(fliplr(fine_bcnts))<= (N(2)-pnt_cnt));
    pnt_cnt = pnt_cnt + sum(fine_bcnts((length(fine_xloc) - closest_ind(end))+1:length(fine_xloc)));
    if abs(N(2) - pnt_cnt) < 1 | ...
        ~sum(fine_bcnts((length(fine_xloc) - closest_ind(end))+1:length(fine_xloc))) | ...
        length(fine_xloc) > closest_ind(end) && N(2) < (pnt_cnt + fine_bcnts((length(fine_xloc) - closest_ind(end)))) %** NEW CONDITION
      max_t = fine_xloc((length(fine_xloc) - closest_ind(end))+1);
    end
    if length(fine_xloc) == 1
      max_t = max(d_vec);
    else
      half_space = (fine_xloc(2)-fine_xloc(1))/2;
    end
    bin_range = fine_xloc(1) + [-half_space,half_space];
    dbn = d_vec(d_vec >= bin_range(1) & d_vec <= bin_range(2));
    nbns = length(dbn);
  end
  
  [l,s] = find(d >= max_t);
  pnt_ind{2} = [l,s];
end
thresh = [min_t, max_t];
dth = d;
dth(dth > max_t) = max_t;
dth(dth < min_t) = min_t;
end