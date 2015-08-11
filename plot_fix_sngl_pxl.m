function plot_fix_sngl_pxl(info_struct,elem_ana_range,varargin);

%% Script:  plot_fix_sngl_pxl(info_struct,elem_ana_range,varargin);
% Description:  Generates a histogram plot of the magnitudes of errors for
% all isotope data denoising was performed on.
% Example:
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% info_struct:  data structure which is the second output from
% fix_sngl_pxl.m (v4)
% varargin - 'PropertyName','PropertyValue'
%   'hdim': numeric indicating whether to plot histograms in 2-D (one
%   histogram plot/isotope) or 3D (one histogram plot for all isotopes)
%   [DEFUALT = 2]
%   'nbins': number of bins in the histogram. The interval 0-1 is evenly
%   divided into nbins. [DEFAULT = 50]
%   'style': bar graph plotting style. See additional documentation and
%   options in help bar.m [DEFAULT = 'stacked']
%   'type': string indicating type of plot=> 'line', or 'bar' [DEFAULT =
%   'line']
%   'id_str': string indicating how user would like to identify data
%   'normalize': binary indicating whether data should be displayed as
%   # of occurances (0) or normalized to percentage of total data (1)
%   [DEFAULT = 0]
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  20 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1

fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(info_struct,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('hdim',PropertyNames)
  hdim = PropertyVal{strmatch('hdim',PropertyNames)};;
else
  hdim = 2;
end

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};;
else
  nbins = 50;
end

if strmatch('style',PropertyNames)
  bar_style = PropertyVal{strmatch('style',PropertyNames)};;
else
  bar_style = 'stacked';
end

if strmatch('type',PropertyNames)
  type = PropertyVal{strmatch('type',PropertyNames)};;
else
  type = 'line';
end

if strmatch('id_str',PropertyNames)
  id_str = sprintf('Sample %s',PropertyVal{strmatch('id_str',PropertyNames)});
else
  id_str = '';
end

if strmatch('norm',PropertyNames)
  nrm =PropertyVal{strmatch('norm',PropertyNames)};
else
  nrm = 0;
end

if strmatch('error',PropertyNames)
  if PropertyVal{strmatch('error',PropertyNames)};
    for f = elem_ana_range
      d = getfield(info_struct,fldmn{f});
      plot_heatmap(d.dfrac,sprintf('%s Magnitude Pxl Change',fldnm{f}),'plot','lin');
    end
  end
end

title_str = sprintf('%s Noise Points Fractional Difference',id_str);

%% Plots Image error as a bar graph (One 2D image / isotope or One 3D image
% containing information for all isotopes
switch type
  case 'bar'
    b = [0:1/nbins:49/nbins] + 1/(2*nbins); % vector containing bin centers for 50 evenly spaced bins between 0 and 1
    %if hdim == 2
    for f = elem_ana_range
      iso_dat = getfield(info_struct,fldnm{f});
      dfrac = getfield(iso_dat,'dfrac');
      dfracvec = dfrac((iso_dat.fxed_all(2,:)-1).*size(dfrac,1)+iso_dat.fxed_all(1,:));
      d10 = dfrac((iso_dat.fxed_10(2,:)-1).*size(dfrac,1)+iso_dat.fxed_10(1,:));
      d100 = dfrac((iso_dat.fxed_100(2,:)-1).*size(dfrac,1)+iso_dat.fxed_100(1,:));
      d110 = dfrac((iso_dat.fxed_110(2,:)-1).*size(dfrac,1)+iso_dat.fxed_110(1,:));
      x10 = hist(d10,b);
      x110 = hist(d110,b);
      x100 = hist(d100,b);
      figure('color','w');
      bar(b,[x10',x110',x100'],bar_style);
      hold on;
      legend({'high-low','high-low-low','high-high-low'},'fontsize',fntsz);
      title(sprintf('%s %s',fldnm{f},title_str),'fontsize',fntsz)
      set(gca,'fontsize',14,'OuterPosition',[0,0.1,1,0.9])
      xlabel('$$\frac{d_{old} - d_{new}}{d_{old}}$$','fontsize',fntsz,'interpreter','latex')
      ylabel('Number of points changed','fontsize',fntsz);
      axis tight
    end
    %>>>>>>>>>>>>>>>>DEBUGGING AS OF 9/21 - See makebars.m to understand
    %requried data structure for a 3D stacked bar plot
    %   else
    %     x10 = zeros(nbins,length(elem_ana_range));
    %     x110 = zeros(nbins,length(elem_ana_range));
    %     x100 = zeros(nbins,length(elem_ana_range));
    %     z_str = cell(0,length(elem_ana_range));
    %     i = 0;
    %     for f = elem_ana_range
    %       %Stacked 3-D bar graphs
    %       %http://www.mathworks.cn/matlabcentral/newsreader/view_thread/258015
    %       i = i+1;
    %       iso_dat = getfield(info_struct,fldnm{f});
    %       dfrac = getfield(iso_dat,'dfrac');
    %       dfracvec = dfrac((iso_dat.fxed_all(2,:)-1).*size(dfrac,1)+iso_dat.fxed_all(1,:));
    %       d10 = dfrac((iso_dat.fxed_10(2,:)-1).*size(dfrac,1)+iso_dat.fxed_10(1,:));
    %       d100 = dfrac((iso_dat.fxed_100(2,:)-1).*size(dfrac,1)+iso_dat.fxed_100(1,:));
    %       d110 = dfrac((iso_dat.fxed_110(2,:)-1).*size(dfrac,1)+iso_dat.fxed_110(1,:));
    %       x10(:,i) = hist(d10,b);
    %       x110(:,i) = hist(d110,b);
    %       x100(:,i) = hist(d100,b);
    %       z_str{i} = fldnm{f};
    %     end
    %     figure('color','w');
    %     bh = bar3(b,x110);
    %
    %     for i=1:length(bh)
    %       zz = get(bh(i),'Zdata');
    %       k = 1;
    %       % Bars are defined by 6 faces(?), adding values from data2 will
    %       % shift the bars upwards accordingly, I'm sure this could be made
    %       % better!
    %       for j = 0:6:6*(length(bh)-1)
    %         zz(j+1:j+6,:)=zz(j+1:j+6,:)+x100(k,i)+x10(k,i);
    %         k=k+1;
    %       end
    %       %---> How MATLAB does stacked bars in makebars.m
    %       %if plottype==1 && m>1, % Stacked
    %       %    z = cumsum(z.').';
    %       %    zz = zz + [zeros(nn,4) z(ones(6,1)*(1:n),ones(4,1)*(1:m-1))];
    %       %end
    %       % Reset Zdata in chart
    %       set(bh(i),'Zdata',zz);
    %     end
    %     set(bh,'FaceColor',[1 0 0]);
    %     % Apply hold so that data2 can be plotted
    %     hold on;
    %     % Plot data2
    %     bh=bar3(data2);
    %     % Set face color to blue
    %     set(bh,'FaceColor',[0 0 1]);
    %     hold off;
    %
    %     hold on;
    %     legend({'high-low','high-low-low','high-high-low'},'fontsize',fntsz);
    %     title(sprintf('%s 041211\\Abrain frac diff',fldnm{f}),'fontsize',fntsz)
    %     xlabel('$$\frac{d_{old} - d_{new}}{d_{old}}$$','fontsize',fntsz,'interpreter','latex')
    %     ylabel('Number of points changed','fontsize',fntsz);
    %     set(gca,'fontsize',14)
    %     axis tight
    %   end
  case 'line'
    b = [0:1/nbins:49/nbins] + 1/(2*nbins);
    cnts = zeros(nbins,length(elem_ana_range));
    %x10 = zeros(nbins,length(elem_ana_range));
    %x110 = zeros(nbins,length(elem_ana_range));
    %x100 = zeros(nbins,length(elem_ana_range));
    z_str = cell(0,length(elem_ana_range));
    i = 0;
    for f = elem_ana_range
      if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
        i = i+1;
        iso_dat = getfield(info_struct,fldnm{f});
        dfrac = getfield(iso_dat,'dfrac');
        if ~isempty(iso_dat.fxed_all)
          dfracvec = dfrac((iso_dat.fxed_all(2,:)-1).*size(dfrac,1)+iso_dat.fxed_all(1,:));
          cnts(:,i) = hist(dfracvec,b);
          y_str = 'Number of points changed';
          if nrm
            cnts(:,i) = (cnts(:,i)./ size(iso_dat.fxed_all,2)).*100;
            y_str ='% points';
          end
          %d10 = dfrac((iso_dat.fxed_10(2,:)-1).*size(dfrac,1)+iso_dat.fxed_10(1,:));
          %d100 = dfrac((iso_dat.fxed_100(2,:)-1).*size(dfrac,1)+iso_dat.fxed_100(1,:));
          %d110 = dfrac((iso_dat.fxed_110(2,:)-1).*size(dfrac,1)+iso_dat.fxed_110(1,:));
          %x10(:,i) = hist(d10,b);
          %x110(:,i) = hist(d110,b);
          %x100(:,i) = hist(d100,b);
          z_str{i} = sprintf('%s (%1.2f%%), d_f %1.2e',fldnm{f},iso_dat.perc_pnts,iso_dat.dfracsum);
        else
          i = i - 1;
        end
      end
    end
    cnts = cnts(:,1:i);
    z_str = z_str(1:i);
    figure('color','w');
    hold on;
    colr = colormap(jet(i));
    for j = 1:i
      plot(b,cnts(:,j),'linewidth',2,'color',colr(j,:));
    end
    legend(z_str,'fontsize',fntsz);
    title(title_str,'fontsize',fntsz)
    set(gca,'fontsize',14,'OuterPosition',[0,0.1,1,0.9])
    xlabel('$\frac{d_{old} - d_{new}}{d_{old}}$','fontsize',fntsz+2,'interpreter','latex')
    ylabel(y_str,'fontsize',fntsz);
    axis tight
end

dockf on all
end

