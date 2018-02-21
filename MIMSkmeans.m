function [] = MIMSkmeans(I)
%% Script: [thresh,pnt_ind] = auto_thresh(d,modestr,nbins)
% Description: Backward search to find the intensity data that makes up the highest x-% of values
% Example:
% INPUTS ----------------------------------------------------------------
% d:  l x s matrix containing intensity values for the given element (l = #
% lines, s = # samples)
% modestr:  1 x 2 matrix where [%low_intensities, %high_intensities] are to be
% eliminated to define thresholds min and max thresholds. This can also be
% indicated as a string 'auto' where the default is used
%   [DEFAULT = bottom 0.02% and top 0.5%]
%     INPUT OPTIONS FOR x:
%       • 'auto': [DEFUALT] thresholds based on a coarse histogram
%       • 'smalliqr': defines image threshold boundaries by defining
%       outliers using the looser statistical definition:
%       data <= mean - 1.5(IQR) and >= mean + 1.5(IQR)  are outside the
%       display range (where IQR = inter quartile range)
%       • 'largeiqr': defines image threshold boundaries by defining
%       outliers using the more strict statistical definition:
%       data <= mean - 3(IQR) and >= mean + 3(IQR)  are outside the
%       display range (where IQR = inter quartile range)
%       • 'cutinfreq': cut at infrequently occurring values. looks at the
%       unique values in the data and the frequency of occurrence for each
%       of the unique values. If a block at the beginning or the end of the
%       histogram occur less than 2 times for an extended block, the data
%       will be thresholded at these levels
%       • 
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
%   1 Nov  2015   AG Balderrama     amanda.gaudreau@gmail.com     4.1
%     Added interquartile range criteria for outlier thresholds and
%     definition of cutoff values
%  24 Feb  2016   AG Balderrama     amanda.gaudreau@gmail.com     4.2
%     Adding an option to cut low counts

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
  binmethod = 0;
end

if ~exist('modestr') || isempty(modestr)
  modestr = 'auto';
end

if ~isempty(strmatch('auto',modestr)) || ~isempty( strmatch('block',modestr))
  x = [2e-2,5e-1];
elseif ~isstr(modestr)
  x = modestr;
end

d_vec = d(:);
if ~isempty(strmatch('smalliqr',modestr))
  q1 = quantile(d(:),0.25);
  q3 = quantile(d(:),0.75);
  iqrd = q3 - q1;
  if iqrd == 0
    iqrd = nanmean(unique(d(:)));
  end
  min_t = max([0,q1-1.5*iqrd]);
  max_t = q3 + 1.5 * iqrd;
elseif ~isempty(strmatch('largeiqr',modestr))
  q1 = quantile(d(:),0.25);
  q3 = quantile(d(:),0.75);
  iqrd = q3 - q1;
  if iqrd == 0
    iqrd = std(d(:));
  end
  min_t = max([0,q1-3*iqrd]);
  max_t = q3 + 3 * iqrd;
elseif ~isempty(strmatch('cut',modestr))
  [ud,~,udin] = unique(d(~isnan(d)));
  hud = hist(udin,1:length(ud));
  hubin = medfilt1( double(hud <= 2),7);
  hbwcc = bwconncomp(hubin);
  if ismember(1,hbwcc.PixelIdxList{1})
    min_t = ud(hbwcc.PixelIdxList{1}(end)+1);
  else
    min_t = ud(1);
  end
  
  if ismember(max(udin),hbwcc.PixelIdxList{end})
    max_t = ud(hbwcc.PixelIdxList{end}(1)-1);
  else
    max_t = ud(end);
  end
elseif ~isempty(strmatch('block',modestr))
  thresh = quantile(d_vec,[x(1)*1e-2,1-x(2)*1e-2]);
  %-------------> New code for 'auto' method added 24 Feb 2016
  [ud,~,udin] = unique(d(~isnan(d)));
  hud = hist(udin,1:length(ud));
  %figure; subplot(1,2,1); plot(hud); subplot(1,2,2); plot(hubin,'.');
  hubin = medfilt1(double(hud <= min([mean(hud),2])),3);
  hbwcc = bwconncomp(hubin);
  plist = cellfun(@(c) sum(ismember(c,find(ud >= thresh(2)))),hbwcc.PixelIdxList);
  pind = find(plist > 0);
  if ~isempty(pind) && ud(hbwcc.PixelIdxList{pind(1)}(1)) < thresh(2) && ...
      ud(hbwcc.PixelIdxList{pind(1)}(1)) > thresh(1)
    thresh(2) = ud(hbwcc.PixelIdxList{pind(1)}(1));
  end
  %<-------------
  min_t = thresh(1); max_t = thresh(2);
elseif binmethod == 0 %% CURRENT DEFAULT MODE - this code is execited when modestr = 'auto' or is not indicated
  thresh = quantile(d_vec,[x(1)*1e-2,1-x(2)*1e-2]);
  min_t = thresh(1); max_t = thresh(2);
else
  error = 0;  % Can be set to a non-zero value
  npnts = size(d,1)*size(d,2);
  N = modestr.*npnts/100;   % Translates percentage of data in number of points, N
  
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
end
[l,s] = find(d <= min_t);
pnt_ind{1} = [l,s];
[l,s] = find(d >= max_t);
pnt_ind{2} = [l,s];
thresh = [min_t, max_t];
dth = d;
dth(dth > max_t) = max_t;
dth(dth < min_t) = min_t;
end