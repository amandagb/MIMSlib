function bkgnd_tissue_ana(data,elem_ana_range,varargin);
%% Script:
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'color_res':  numeric specifying color resolution of the heatmap
%   (maximum 256) [DEFAULT = 256]
%   'thresh_data':  1 x 2 matrix with minimum threshold value as first
%   element and maximum threshold value as second
%   'lines': matrix containing line numbers user would like to plot
%   [DEFAULT: ALL LINES]
%   'samples': matrix containing samples (S) user would like to plot [DEFAULT:
%   ALL SAMPLES]
%   'plot_type': string - 'log','linear','both' [DEFUALT: BOTH]
%   'elim_blanks': 1 or 0 - eliminates any lines which have NAN data
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  30 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('norm',PropertyNames)
  norm_str = PropertyVal{strmatch('norm',PropertyNames)};
else
  norm_str = 'lpfmed';
end

if strmatch('nbin',PropertyNames)
  nbin = PropertyVal{strmatch('nbin',PropertyNames)};
else
  nbin = 1e2;
end

[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
[dnorm,vco] = drift_correct(data,elem_ana_range,'plot',0,'stat',0,'freq',0,varargin{:});
hco = find_horiz_co(data,elem_ana_range,'dnorm',dnorm,'vco',vco,'plot',0);
for f = elem_ana_range
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    if isstruct(data)
      d = getfield(dnorm,fldnm{f});
      v = getfield(vco,fldnm{f});
      h = getfield(hco,fldnm{f});
    end
    upco = 8;%h.up_co;
    lowco = 125;%h.low_co;
    lco = 82;%v.l_co;
    rco = 402;%v.r_co;
    la0 = 505:525;%v.la0ind;
    bk = 1:la0(1)-1;
    
    
    top = d(1:upco,bk);
    bottom = d(lowco:size(d,1),bk);
    left = d(upco+1:lowco-1,1:lco);
    right = d(upco+1:lowco-1,rco:bk(end));
    dbkgnd = [top(:);bottom(:);left(:);right(:)];
    Nbk = length(dbkgnd);
    
    dla0 = d(:,la0);
    dla0 = dla0(:);
    Nla0 = length(dla0);
    
    dtissue = d(upco+10:lowco-10,lco+10:rco-10);
    dtissue = dtissue(:);
    Ntissue = length(dtissue);
    
    figure('color','w');

    [t,ind] = auto_thresh(dtissue,[0.1,0.5],[]);
    dtissue(dtissue > t(2)) = t(2);
    dtissue(dtissue < t(1)) = t(1);
    [cnttissue,xtissue] = hist(dtissue,nbin);
    bar(xtissue,cnttissue./max(cnttissue),'barwidth',1,'facecolor','r','displayname','Tissue')
    
    hold on;
    [t,ind] = auto_thresh(dbkgnd,[0.1,0.5],[]);
    dbkgnd(dbkgnd > t(2)) = t(2);
    dbkgnd(dbkgnd < t(1)) = t(1);
    [cntbk,xbk] = hist(dbkgnd,nbin);
    bar(xbk,cntbk./max(cntbk),'barwidth',1,'facecolor','b','displayname','Background')
    
    %[t,ind] = auto_thresh(dla0,[0.1,0.5],[]);
    %dla0(dla0 > t(2)) = t(2);
    %dla0(dla0 < t(1)) = t(1);
    [cntla0,xla0] = hist(dla0,50);
    bar(xla0,cntla0./max(cntla0),'barwidth',1,'facecolor','w','displayname','Laser Off')
    
    legend('Background','Laser Off','Tissue');    
    title(sprintf('Normalized histograms for %s',fldnm{f}),'fontsize',14);
    xlabel('Intensity (cps)','fontsize',14);
    ylabel('Normalized # of points','fontsize',14);
    dockf on all
  end
end
