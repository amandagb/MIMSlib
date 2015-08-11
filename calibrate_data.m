function lnparam = calibrate_data(data, elem_ana_range,varargin)

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
%  3  Oct  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

% Parameters that are specific to 'Standards\032211' data
ngroups = 9;
lines_per_conc = 3;
nconc = 7;
relevant_col = 100:175;
conc = [0, 3.2, 6.5, 12.9, 19.4, 25.8, 32.3];

% Plotting variables
colors = cool(ngroups);
line_marker = {'o','x','s'};
fntsz = 14;

lnparam = []
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 1e3;
end

[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
for f = elem_ana_range
  if isstruct(data)
    d = getfield(data,fldnm{f});
  else
    d = data;
  end
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    ind = 1:(ngroups*lines_per_conc*nconc);
    %start_ind = 1:lines_per_conc:(ngroups*lines_per_conc*nconc);
    %end_ind = start_ind(2:end)-1;
    m_conc_data = zeros(ngroups*lines_per_conc,nconc);
    med_conc_data = zeros(ngroups*lines_per_conc,nconc);
    for c = 1:nconc
      v = repmat([zeros(1,lines_per_conc*(c-1)),ones(1,lines_per_conc),zeros(1,(lines_per_conc*(nconc - c)))],1,ngroups);
      conc_lines = ind(v==1);
      conc_data = d(conc_lines,relevant_col);
      m_conc_data(:,c) = mean(conc_data,2);
      med_conc_data(:,c) = median(conc_data,2);
    end
    
    figure('color','w','name',fldnm{f});
    hold on;
    for g = 1:ngroups
      for i = 1:3
        if strmatch(lower(fldnm{f}),'au')
          plot(g.*ones(length(conc),1),m_conc_data((g-1)*lines_per_conc+i,:),'marker',line_marker{i},'color',colors(g,:),...
%           plot([1:5,7,6],m_conc_data((g-1)*lines_per_conc+i,:),'marker',line_marker{i},'color',colors(g,:),...
          'linestyle','none');
        else
          plot(conc,m_conc_data((g-1)*lines_per_conc+i,:),'marker',line_marker{i},'color',colors(g,:),...
            %           plot([1:5,7,6],m_conc_data((g-1)*lines_per_conc+i,:),'marker',line_marker{i},'color',colors(g,:),...
          'linestyle','none');
        end
      end
    end
    x = repmat(conc,ngroups*lines_per_conc,1);
    x = x(:);
    y = m_conc_data(:);
    %figure; plot(x,y,'.');
    [p,S] = polyfit(x,y,1);
    hold on;
    plot(x,p(1).*x+p(2),'r');
    %plot(conc,mean(m_conc_data,1),'marker','.','color','k','linestyle','none');
    set(gca,'fontsize',fntsz);
    title(sprintf('Calibration curve for %s',fldnm{f}),'fontsize',fntsz);
    xlabel('Concentrations [ppm]','fontsize',fntsz);
    ylabel('Intensities [cps]','fontsize',fntsz);
    lnparam = setfield(lnparam,fldnm{f},p);
  end
end

