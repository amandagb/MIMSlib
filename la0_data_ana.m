function la0_data_ana(data,elem_ana_range,la0ind,varargin);

%% Script:  la0_data_ana(data,elem_ana_range,la0ind,varargin);
% Description:  Generates a histogram plot of the magnitudes of errors for
% all isotope data denoising was performed on.
% Example:  la0_data_ana(dA,[],ila0);
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% elem_ana_range: 1 x <# elements> matrix which corresponds to the field
%                 number in the data structure. User may also specify
%   This variable can also be a cell of strings where the string specifies
%   the isotope the user would like to analyze
% la0ind: 1 x <sample #> row vector indicating the columns (sample #'s) for
% which the LA system was off
% varargin - 'PropertyName','PropertyValue'
%   'hdim': numeric indicating whether to plot histograms in 2-D (one
%   histogram plot/isotope) or 3D (one histogram plot for all isotopes)
%   [DEFUALT = 3]
%   'nbins': number of bins in the histogram. The interval 0-1 is evenly
%   divided into nbins. [DEFAULT = 50]
%   'style': bar graph plotting style. See additional documentation and
%   options in help bar.m [DEFAULT = 'stacked']
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  23 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1

fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('hdim',PropertyNames)
  hdim = PropertyVal{strmatch('hdim',PropertyNames)};
else
  hdim = 2;
end

if strmatch('tit',PropertyNames)
  tit_str = strcat(', ',PropertyVal{strmatch('tit',PropertyNames)});
else
  tit_str = '';
end

if strmatch('data_id',PropertyNames)
  id_str = sprintf('Sample %s',PropertyVal{strmatch('id_str',PropertyNames)});
else
  id_str = '';
end

%% Plots Image error as a bar graph (One 2D image / isotope or One 3D image
% containing information for all isotopes
i = 0;
line_mean = zeros(size(getfield(data,fldnm{1}),1),length(elem_ana_range));
line_med = zeros(size(getfield(data,fldnm{1}),1),length(elem_ana_range));
line_std = zeros(size(getfield(data,fldnm{1}),1),length(elem_ana_range));
%la0dat = zeros(size(getfield(data,fldnm{1}),1),length(elem_ana_range)*length(la0ind));
for f = elem_ana_range
  if isstruct(data)
    d = getfield(data,fldnm{f});
  else
    d = data;
  end
  
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    i = i + 1;
    la0dat = d(:,la0ind);
    line_mean(:,i) = mean(la0dat')';
    line_med(:,i) = median(la0dat')';
    line_std(:,i) = std(la0dat')';
    sample_mean = mean(la0dat(:));
    sample_std = std(la0dat(:));
    
    figure('color','w','name',fldnm{f})
    %plot([line_mean(:,i),line_med(:,i)]);
    plot(line_med(:,i),'linewidth',2);
    hold on;
    %errorbar(line_mean(:,i),line_std(:,i),'r');
    plot(line_mean(:,i),'r','linewidth',2);
    plot([line_mean(:,i) + line_std(:,i),line_mean(:,i) - line_std(:,i)],':r');
    %ylim([sample_mean - 2*sample_std,sample_mean + 2*sample_std])
    axis tight
    plot(xlim,sample_mean.*[1,1],'c');
    plot(xlim,(sample_mean + sample_std).*[1,1],':c');
    plot(xlim,(sample_mean -sample_std).*[1,1],':c');
    title({[sprintf('Laser off statistics for %s',fldnm{f})],[sprintf('%s %s',tit_str, id_str)]}...
      ,'fontsize',fntsz);
    xlabel('Line #','fontsize',fntsz);
    ylabel('Mean Intensity (cps)','fontsize',fntsz);
    legend({'median','mean','standard deviation'},'fontsize',fntsz,'location','sw');
    set(gca,'fontsize',fntsz);
  end
end

% figure('color','w')
% plot(line_mean,'linestyle','-','marker','.');
% hold on;
% plot(line_med,'linestyle','-','marker','*');
% %plot(line_std,'linestyle','-','marker','*');
dockf on all

end

