function out = plot_hist(data,elem_ana_range,varargin)

%% Script:  plot_hist(data,elem_ana_range,'nbins',#,'thresh',[low,high],'perc',[],'norm',[],'lines',[low:high],'samples',[low,high],'xzoom',[low,high],'yzoom',[low,high])
% Description: Plots a histogram of the data
% Example:
% plot_hist(data,'zn','thresh',[0,4e6],'samples',200:300,'lines',60:120,'nbins',1e2);
% Required Function: interp_data_elemrng,
%
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrice for each
%       isotype analyzed
% elem_ana_range: 1 x <# elements> matrix which corresponds to the field
%                 number in the data structure. User may also specify
%   This variable can also be a cell of strings where the string specifies
%   the isotope the user would like to analyze
% varargin - 'PropertyName','PropertyValue'
%   'nbins':  numeric specifying the number of histogram bins
%   'thresh_data':  1 x 2 matrix with minimum threshold value as first
%       element and maximum threshold value as second
%   'percent':  []. Takes <# lines>*<#samples> and divides ever hist bin
%       count by that number
%   'normalize':  []. Takes the hist bin with the maximum count and divides
%       every count by that number
%   'lines': matrix containing line numbers user would like to plot
%       [DEFAULT: ALL LINES]
%   'samples': matrix containing samples (S) user would like to plot [DEFAULT:
%       ALL SAMPLES]
%   'xzoom': 1 x 2 matrix with lower and upper x-limit
%   'yzoom': 1 x 2 matrix with lower and upper y-limit
%
% OUTPUTS ---------------------------------------------------------------
% out:  structure output. fields can be added as user find output data which
% are desirable/necessary
%
%  Date           Author            E-mail                      Version
%  21 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

%% Initialize user defined variables
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 1e3;
end

int_range = [];
if strmatch('int_range',PropertyNames)
  int_range = PropertyVal{strmatch('int_range',PropertyNames)};
end

fi = 0;
out = [];
for f = elem_ana_range
  h = zeros(length(elem_ana_range),1);
  fi = fi+1;
  if ~exist('d','var')
    d = getfield(data,fldnm{f});
  end
  if strmatch(fldnm{f},'Time')
    fstr = sprintf('%s Difference Histogram', fldnm{f});
    h(fi)=figure('Name',fstr,'color','w','position',[360,502,560,420]);
    subplot(1,2,1)
    tdiff = diff(d);
    hist(tdiff(:),nbins)
    xlabel('Time Difference [sec]');
    ylabel('# of occurrences');
    title('Time difference between same sample points on different lines');
    
    subplot(1,2,2)
    tdiff = diff(d,1,2);
    hist(tdiff(:),nbins)
    xlabel('Time Difference [sec]');
    ylabel('# of occurrences');
    title('Time difference between next sample point on same lines');
  elseif isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    elem_data = [];
    elem_data = setfield(elem_data,'data',d);
    elem_str = strcat(fldnm{f},' Histogram');
    
    thresh_str = '';
    if strmatch('thresh',PropertyNames)
      thresh_data = PropertyVal{strmatch('thresh',PropertyNames)};
      thresh_d = d(:);
      if ~isempty(strmatch('auto',thresh_data)) || (sum(thresh_data - 10) > 0)
        thresh_data = auto_thresh(d,thresh_data,[]);
      end
      ind_gtT = find(thresh_d > max(thresh_data));
      ind_ltT = find(thresh_d < min(thresh_data));
      thresh_d(ind_gtT) = max(thresh_data);
      thresh_d(ind_ltT) = min(thresh_data);
      d = reshape(thresh_d,size(d));
      thresh_str = sprintf('THRESHOLDS: lower = %g; upper =  %g; ',min(thresh_data),max(thresh_data));
    end
    
    line_str = '';
    if strmatch('lines',PropertyNames)
      l = PropertyVal{strmatch('lines',PropertyNames)};
      d = d(l,:);
      line_str = strcat('LINES: ',num2str(min(l)),'-',num2str(max(l)));
    end
    
    sample_str = '';
    if strmatch('samples',PropertyNames)
      s = PropertyVal{strmatch('samples',PropertyNames)};
      d = d(:,s);
      sample_str = strcat('SAMPLES: ',num2str(min(s)),'-',num2str(max(s)));
    end
    
    elim_str = '';
    if strmatch('elim',PropertyNames)
      l = PropertyVal{strmatch('elim',PropertyNames)};
      [nan_ln,v] = find(isnan(d)==1);
      nan_ln = unique(nan_ln);
      d(nan_ln,:) = [];
      elim_str = strcat(num2str(length(nan_ln)),' lines eliminated');
    end
    
    fstr = sprintf('%s Histogram', fldnm{f});
    h(fi) = figure('Name',fstr,'color','w','position',[360,502,560,420]);
    [bin_cnts,x_loc] = hist(d(:),nbins);
    bin_sz = x_loc(2) - x_loc(1);
    bar(x_loc,bin_cnts,1);
    xlabel('Intensity [cps]','fontsize',12);
    ylabel('# of occurrences','fontsize',12);
    
    if strmatch('perc',PropertyNames)
      [bin_cnts,x_loc] = hist(d(:),nbins);
      tot_pnts = sum(bin_cnts);
      bar(x_loc,(bin_cnts./tot_pnts).*100,1);
      xlabel('Intensity [cps]','fontsize',12);
      ylabel('% of occurences','fontsize',12);
    end
    
    if strmatch('norm',PropertyNames)
      [bin_cnts,x_loc] = hist(d(:),nbins);
      bar(x_loc,bin_cnts./max(bin_cnts),1);
      xlabel('Intensity [cps]','fontsize',12);
      ylabel('Normalized # of occurences','fontsize',12);
    end
    
    title({strcat(elem_str),...
      [strcat(thresh_str,sprintf(' BIN SIZE: %g',bin_sz))],...
      [strcat(line_str,sample_str)]},'fontsize',12);
    if strmatch('xzoom',PropertyNames)
      s = PropertyVal{strmatch('xzoom',PropertyNames)};
      xlim(s);
    end
    
    if strmatch('yzoom',PropertyNames)
      s = PropertyVal{strmatch('yzoom',PropertyNames)};
      ylim(s);
    end
    
    elem_data = setfield(elem_data,'data_plotted',d);
    elem_data = setfield(elem_data,'bin_cnts',bin_cnts);
    elem_data = setfield(elem_data,'x_loc',x_loc);
    out = setfield(out,fldnm{f},elem_data);
    clear elem_data d
  end
end
dockf on all
end