function varargout = plot_heatmap(data,elem_ana_range,varargin)

%% Script: [varargout{figure handle, plotted data, x coords, y coords}]=
%   plot_heatmap(data,elem_ana_range,varargin)
% Description:
% Example:  plot_heatmap(d,{'cu','zn'},'thresh','auto')
%           plot_heatmap(d,[3,6])
%           [h,d,x,y] = plot_heatmap(d,[]);
% [h,d,x,y] = plot_heatmap(d,[],'color_res',256,'thresh',0,...
%   'lines',[1:size(data,1)], 'samples', [1:size(data,2)],...
%   'plot_type','lin','elim_blanks',0,'scale',0,'overlay',0,...
%   'handle',0,'type','i','showtitle',1,'showcbar,1,'title_txt','title text',...
%   'channel_order',{'rgbymcw'},'save_heat',0,'save_str','file_name');
% Required Function: interp_data_elemrng, custom_colormaps
%
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% elem_ana_range: 1 x <# elements> matrix which corresponds to the field
%                 number in the data structure. User may also specify
%   This variable can also be a cell of strings where the string specifies
%   the isotope the user would like to analyze
% varargin - 'PropertyName','PropertyValue'
%   • 'color_res':  numeric specifying color resolution of the heatmap
%   (maximum 256) [DEFAULT = 256]
%   • 'thresh_data':  1 x 2 matrix with minimum threshold value as first
%   element and maximum threshold value as second. Can also use string
%   'auto' to cut off bottom 0.02% and top 0.5% of data points
%   • 'lines': matrix containing line numbers user would like to plot
%   [DEFAULT: ALL LINES]
%   • 'samples': matrix containing samples (S) user would like to plot [DEFAULT:
%   ALL SAMPLES]
%   • 'plot_type': string - 'log','linear','both' [DEFUALT: 'lin']
%   • 'elim_blanks': 1 or 0 - eliminates any lines which have NAN data
%   • 'scale': 1 or 0 - uses spatial scaling rather than line # and sample #
%   (1) or not (0) [DEFAULT = 0]
%   • 'spot_size': integer indicating spot size in microns
%   • 'l2l_dist': integer indicating center to center spacing of lines, in microns
%   • 'sampling_time': integer indicating time between successive samples, in
%   seconds
%   • 'scan_speed': integer indicating scan speed in microns/sec
%   • 'xlparam': logical indicating whether to use 'CBM Data Log v3.xls'
%   file or not [DEFAULT = 0]
%   • 'rottag': rotation tag (as defined in rotMIMS.m file) [DEFAULT = 0]
%         -2    Original data aquired VRDL, flip columns up/down (or flip LR
%               then rotate twice)
%         -1    Original data aquired RVLD, rotate 90 deg clockwise
%          0    Original data aquired LDRV, rotate 90 deg ccw
%          1    Original data aquired LDRV, rotate 90 deg ccw
%          2    Original data aquired VLDR, rotate 180 deg
%   • 'rotreq': logical indicating whether rotation is required (1) or not
%   (0) [DEFAULT = 1]
%   • 'overlay': binary indicating whether to overlay input isotopes (1) or
%   not (0) [DEFAULT = 0]
%   • 'handle': figure handle user would like to generate the images on. (0)
%   to generate a new figure. [DEFAULT = 0]
%   • 'type': string indicating whether data is given in insensity (counts
%   per second) -- 'i'-- or absolute concentration (parts per million)
%   --'c'-- [DEFAULT = 'i']
%   • 'showtitle': binary indicating whether to display the title on the
%   final image (1) or not (0) [DEFAULT = 1]. In addition, strings 'd' and
%   'h' can be used to display the plotted details only or the header only,
%   respectively.
%   • 'title_txt': string indicating the preferred text. Element name will
%   also be included
%   • 'data_mod': string which is a data modifier. For instance, when
%   plotting denoised data, this string could be indicated as 'denoised'
%   and will be added to the figure name and title strings.
%   • 'channel_order': vector the same length as elem_ana_range indicating
%   the order of colors (1 = R, 2 = G, 3 = B, 4 = Y, 5 = M, 6 = C, 7 = W, 8 = jet, 9 = fire, 10 = lum).
%   Used when overlay value == 1.
%   • 'save_heat': binary indicating whether to save heatmaps (1) or not (0)
%   [DEFAULT = 0]
%   • 'save_str': any additional string that should be used to describe
%   figures
%   • 'nomarks': binary indicating whether to remove all marks from the plot
%   (1) or not (0) [DEFAULT = 0]
%   • 'showlegend': binary indicating whether to show the legend (1) or
%   not (0) [DEFAULT = 0] -- should be indicated to override nomarks (can
%   be used with 'nomarks',1,'showleg',1 in order to still show the legend)
%   • 'showcbar': binary indicating whether to show the colorbar (1) or
%   not (0) [DEFAULT = 1]
%   • 'cbarloc': string indicating the colorbar location [DEFAULT = 'westoutside']
%   • 'dockall': binary indicationg whether to dock all figures (1) or not
%   (0) [DEFAULT = 0]
%   • 'figpos': 4 x 1 vector indicating the figure position [DEFAULT = []]
%   • 'normalizeColor': binary indicating whether to normalize color
%   (cps/area[um^2]) (1) or not (0) [DEFAULT = 0] This property is only
%   used if both scale is true. If true, type string will be changed
%   reflect normalization units (cps/area)
%   • 'addlabels': cell of strings <# lines> x 1 indicating the label
%   associated with each line. IF value has more than 1 column, the input
%   data will assumed to the hdrtxt structure
%   • 'useMask': binary matrix the same size as the data which indiates
%   which pixels the user would like to plot. Pixels to be plotted should
%   be indicated by a 1/true and others by 0 [DEFAULT = true(M,N)]
%
% OUTPUTS ---------------------------------------------------------------
% varargout
%   {1}: Figure handle
%   {2}: Plotted data structure (if overlay==0) or matrix (if overlay==1)
%   {3}: Scaling/number of x-axis
%   {4}: Scaling/number of y-axis
%
%  Date           Author            E-mail                      Version
%   7 Apr  2015   AGBalderrama   amandagbalderrama@gmail.com      8
%     Allows user to specify map labels to plot or to add map label names
%     to the image

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

%--------------------------------------------------
%% Initialize variables that deal with scaling the image (converting pixels --> cm)
%--------------------------------------------------
if strmatch('rottag',PropertyNames)
  rotTag = PropertyVal{strmatch('rottag',PropertyNames)};
else
  rotTag = 0; %microns
end

if strmatch('rotreq',PropertyNames)
  rotreq = PropertyVal{strmatch('rotreq',PropertyNames)};
else
  rotreq = 1; %microns
end

if strmatch('xlparam',PropertyNames)
  readxl = PropertyVal{strmatch('xlparam',PropertyNames)};
else
  readxl = 0; %microns
end

if strmatch('spot',PropertyNames)
  spot_size = PropertyVal{strmatch('spot',PropertyNames)};
else
  spot_size = 20; %microns
end

if strmatch('l2l',PropertyNames)
  l2l_dist = PropertyVal{strmatch('l2l',PropertyNames)};
else
  l2l_dist = 20; %center to center line spacing, microns
end

if isstruct(data) && strmatch('Time',fieldnames(data))
  if strmatch('samp',PropertyNames)
    t_samp = PropertyVal{strmatch('samp',PropertyNames)};
  else
    t_samp = [];
  end
  
  if isempty(t_samp)
    t_samp = nanmedian(nanmedian(diff(data.Time')));
    if t_samp == 0
      t_samp = nanmedian(nanmedian(diff(data.Time)));
    end
  end
elseif strmatch('samp',PropertyNames)
  t_samp = PropertyVal{strmatch('samp',PropertyNames)};
else
  t_samp = 0.5; % Time between successive samples
end

if strmatch('scan',PropertyNames)
  scan_speed = PropertyVal{strmatch('scan',PropertyNames)};
else
  %scan_speed = spot_size/t_samp;
  scan_speed = spot_size*2;
end

if strmatch('scale',PropertyNames)
  sc = PropertyVal{strmatch('scale',PropertyNames)};
else
  sc = 0;
end

if sc && readxl
  [spot_size, l2l_dist, scan_speed, E, vgas, t_samp] = getphysicalparams(data.dataset);
end

%--------------------------------------------------
%% Figure customization properties
%--------------------------------------------------
if strmatch('showaqarrow',PropertyNames)
  showaqarrow = PropertyVal{strmatch('showaqarrow',PropertyNames)};
else
  showaqarrow = 0;
end

if strmatch('showtitle',PropertyNames)
  showtitle = PropertyVal{strmatch('showtitle',PropertyNames)};
else
  showtitle = 1;
end

if strmatch('title_txt',PropertyNames)
  title_text = PropertyVal{strmatch('title_txt',PropertyNames)};
end

if strmatch('data_mod',PropertyNames)
  dmod = PropertyVal{strmatch('data_mod',PropertyNames)};
else
  dmod = [];
end

if strmatch('nomarks',PropertyNames)
  nomarks = PropertyVal{strmatch('nomarks',PropertyNames)};
else
  nomarks = 0;
end

if strmatch('dock',PropertyNames)
  dock = PropertyVal{strmatch('dock',PropertyNames)};
else
  dock = 0;
end

if strmatch('figpos',PropertyNames)
  fpos = PropertyVal{strmatch('figpos',PropertyNames)};
else
  fpos = [400,250,560,420];
end

if strmatch('showcbar',PropertyNames)
  cb = PropertyVal{strmatch('showcbar',PropertyNames)};
else
  cb = 1;
end

if strmatch('showleg',PropertyNames)
  legon = PropertyVal{strmatch('showleg',PropertyNames)};
else
  legon = 1;
end

if strmatch('cbarloc',PropertyNames)
  cbl = PropertyVal{strmatch('cbarloc',PropertyNames)};
else
  cbl = 'eastoutside';
end

if strmatch('addlabels',PropertyNames)
  relhtxt = PropertyVal{strmatch('addlabels',PropertyNames)};
else
  relhtxt = [];
end

%--------------------------------------------------
%% Color parameters
%--------------------------------------------------
if strmatch('color',PropertyNames)
  color_res = PropertyVal{strmatch('color',PropertyNames)};
  if color_res > 256
    color_res = 256;
  end
else
  color_res = 256;
end

if strmatch('normalizecolor',PropertyNames)
  colorNorm = PropertyVal{strmatch('normalizecolor',PropertyNames)};
else
  colorNorm = 0;
end

%--------------------------------------------------
%% Overlay parameters
%--------------------------------------------------
if strmatch('over',PropertyNames)
  over = PropertyVal{strmatch('over',PropertyNames)};
  if over
    over_logical = [eye(3);1,1,0;1,0,1;0,1,1;1,1,1]; % R, G, B, Y, M, C, W
    over_col = repmat([0:color_res-1]'./color_res,1,3);%.*ones(color_res,3);
    trans_alpha = 1;
  else
    over = 0;
    trans_alpha = 1;
  end
else
  over = 0;
  trans_alpha = 1;
end

if strmatch('plot',PropertyNames)
  ptype = PropertyVal{strmatch('plot',PropertyNames)};
else
  ptype = 'lin';
end

ptag = cell(0,1);
if strmatch('log',ptype)
  ptag{1} = 'Log';
elseif strmatch('lin',ptype)
  ptag{1} = 'Linear';
elseif strmatch('b',ptype)
  ptag = {'Log','Linear'};
end

if strmatch('save_heat',PropertyNames)
  sv = PropertyVal{strmatch('save_heat',PropertyNames)};
  if strmatch('save_str',PropertyNames)
    svstr = PropertyVal{strmatch('save_str',PropertyNames)};
  else
    svstr = '';
  end
  
else
  sv = 0; %intensity [cps] versus concentration [ppm]
end

% if strmatch('saveheat',PropertyNames)
%   sv = PropertyVal{strmatch('saveheat',PropertyNames)};
%   if strmatch('save_str',PropertyNames)
%     svstr = PropertyVal{strmatch('savestr',PropertyNames)};
%   else
%     svstr = '';
%   end
% else
%   sv = 0; %intensity [cps] versus concentration [ppm]
% end

if strmatch('type',PropertyNames)
  type = PropertyVal{strmatch('type',PropertyNames)};
else
  type = 'i'; %intensity [cps] versus concentration [ppm]
end

if strmatch('c',type)
  type = 'Concentration [ppm]';
else
  type = 'Intensity [cps]';
end

handles = [];
dout = [];

if strmatch('chan',PropertyNames)
  chanord = PropertyVal{strmatch('chan',PropertyNames)};
  chanordstr = chanord;
  chanord = zeros(1,length(chanord));
  chanord(strfind(chanordstr,'r') ) = 1;
  chanord(strfind(chanordstr,'g') ) = 2;
  chanord(strfind(chanordstr,'b') ) = 3;
  chanord(strfind(chanordstr,'y') ) = 4;
  chanord(strfind(chanordstr,'m') ) = 5;
  chanord(strfind(chanordstr,'c') ) = 6;
  chanord(strfind(chanordstr,'w') ) = 7;
  chanord(strfind(chanordstr,'j') ) = 8;
  chanord(strfind(chanordstr,'f') ) = 9;
  chanord(strfind(chanordstr,'l') ) = 10;
elseif over
  chanordstr = 'rgbymcw';
  chanordstr = chanordstr(1:length(elem_ana_range));
  chanord = 1:7;
else
  chanord = 8;
end

n = 0;
for pt = 1:length(ptag)
  for f = elem_ana_range
    n = n + 1;
    if isstruct(data)
      d = getfield(data,fldnm{f});
    else
      d = data;
    end
    
    %% Figure handle
    if strmatch('h',PropertyNames)
      h = PropertyVal{strmatch('h',PropertyNames)};
    else
      h = 0;
    end
    
    line_str = '';
    if strmatch('lines',PropertyNames)
      l = PropertyVal{strmatch('lines',PropertyNames)};
      if ~ischar(d)
        d = d(l,:);
      end
      line_str = strcat('LINES: ',num2str(min(l)),'-',num2str(max(l)),'; ');
    end
    
    sample_str = '';
    if strmatch('samples',PropertyNames)
      s = PropertyVal{strmatch('samples',PropertyNames)};
      d = d(:,s);
      sample_str = strcat('SAMPLES: ',num2str(min(s)),'-',num2str(max(s)),'; ');
    end
    [M0,N0] = size(d);
    
    if strmatch('usemask',PropertyNames)
      maskd = PropertyVal{strmatch('usemask',PropertyNames)};
    else
      maskd = ones(M0,N0);
    end
    d(~maskd) = nan;
    
    if rotreq
      d = rotMIMS(d,rotTag);
    end
    [M,N] = size(d);
    
    if isempty(strmatch(fldnm{f},skip_fields))
      if sc && M0 == M
        y = [0,(M-2)*l2l_dist+spot_size]./1e3;
        if exist('t_samp')
          x = [0,scan_speed.*(N-1)*t_samp]./1e3;
        elseif isstruct(data) && strmatch('Time',fieldnames(data))
          x = [0,scan_speed.*data.Time(1,end)]./1e3;
        end
        xstr = 'Distance [mm]';
        ystr = 'Distance [mm]';
        
        if colorNorm
          type = strcat(type(1:end-1),'/mm2]');
          pxlwidth = (scan_speed*t_samp)./1e3;
          pxlheight = l2l_dist./1e3;
          pxlarea = pxlwidth*pxlheight;
          d = d./pxlarea;
        end
      elseif sc && M0 == N
        x = [0,(M0-2)*l2l_dist+spot_size]./1e3;
        if exist('t_samp')
          y = [0,scan_speed.*(N0-1)*t_samp]./1e3;
        elseif isstruct(data) && strmatch('Time',fieldnames(data))
          y = [0,scan_speed.*data.Time(1,end)]./1e3;
        end
        xstr = 'Distance [mm]';
        ystr = 'Distance [mm]';
        
        if colorNorm
          type = strcat(type(1:end-1),'/mm2]');
          pxlwidth = l2l_dist./1e3;
          pxlheight = (scan_speed*t_samp)./1e3;
          pxlarea = pxlwidth*pxlheight;
          d = d./pxlarea;
        end
      else
        x = 1:size(d,2);
        y = 1:size(d,1);
        xstr = 'Sample #';
        ystr = 'Line #';
      end
      
      thresh_str = '';
      if strmatch('thresh',PropertyNames)
        td = PropertyVal{strmatch('thresh',PropertyNames)};
        thresh_data = td(min(size(td,1),n),:);
        thresh_d = d(:);
        if ~isempty(strmatch('auto',thresh_data)) ...
            || ((sum(thresh_data - 10) < 0) && ~any(thresh_data < 0))...
            || (ischar(td) && isempty(td))
          thresh_data = auto_thresh(d,thresh_data,[]);
        else
          thresh_d(1,1) = thresh_data(1);
          thresh_d(2,1) = thresh_data(2);
        end
        %ind_gtT = find(thresh_d > max(thresh_data));
        %ind_ltT = find(thresh_d < min(thresh_data));
        %thresh_d(ind_gtT) = max(thresh_data);
        %thresh_d(ind_ltT) = min(thresh_data);
        %d = reshape(thresh_d,size(d));
        thresh_str = sprintf('THRESHOLDS: lower = %10.2g; upper =  %10.3g; ',min(thresh_data),max(thresh_data));
        maxcbarlabel = '> ';
        if ~any(d(:) < min(thresh_data)) || ~isempty(strfind(type,'ppm'))
          mincbarlabel = '';
        else
          mincbarlabel = '< ';
        end
      else
        maxcbarlabel = '';
        mincbarlabel = '';
      end
      
      elim_str = '';
      if strmatch('elim',PropertyNames)
        l = PropertyVal{strmatch('elim',PropertyNames)};
        [nan_ln,v] = find(isnan(d)==1);
        nan_ln = unique(nan_ln);
        d(nan_ln,:) = [];
        elim_str = strcat(num2str(length(nan_ln)),' lines eliminated');
      end
      
      data_str = '';
      if isfield(data,'dataset')
        data_str = getfield(data,'dataset');
      end
      
      if over
        fstr = 'Heat'; % figure string
      else
        if isempty(dmod)
          fstr = sprintf('%s Heat', fldnm{f}); % figure string
        else
          fstr = sprintf('%s %s Heat',fldnm{f},dmod); % figure string
        end
      end
      fntsz = 12;%18;
      
      ncbarlabels = 6; % number of color bar labels
      if exist('thresh_data')
        rng = max(thresh_data) - min(thresh_data);
        r1 = min(thresh_data); r2 = max(thresh_data);
        d(1,1) = r1; d(1,2) = r2;
      else
        rng = max(d(:)) - min(d(:));
        r1 = min(d(:)); r2 = max(d(:));
      end
      rngdiv = rng/(ncbarlabels-1);
      pwrof10 = floor(log10(rngdiv));
      if pwrof10 == 1
        rngdiv = round(rngdiv);
      else
        rngdiv_lvl1 = round(rngdiv/10^pwrof10)*10^pwrof10;
        if pwrof10 == 1
          rngdiv_lvl2 = 0;
          pwrof10 = 2;
        else
          rngdiv_lvl2 = round((rngdiv-rngdiv_lvl1)/10^(pwrof10-1))*10^(pwrof10-1);
        end
        rngdiv = rngdiv_lvl1 + rngdiv_lvl2;
      end
      ytick_level = [r1, round((r1)/10^(pwrof10-1))*10^(pwrof10-1) + rngdiv.*[1:ncbarlabels-2], r2];
      yticklab = cell(1,ncbarlabels);
      if floor(log10(max(d(:)))) >= 1 && floor(log10(max(d(:)))) <= 3
        yticklab{1} = sprintf('%s%1.0f',mincbarlabel,ytick_level(1));
        for c = 2:ncbarlabels-1
          yticklab{c} = sprintf('%1.0f',ytick_level(c));
        end
        yticklab{ncbarlabels} = sprintf('%s%1.0f',maxcbarlabel, ytick_level(ncbarlabels));
      elseif floor(log10(max(d(:)))) >= -1 && floor(log10(max(d(:)))) <= 0
        yticklab{1} = sprintf('%s%1.2f',mincbarlabel,ytick_level(1));
        for c = 2:ncbarlabels-1
          yticklab{c} = sprintf('%1.2f',ytick_level(c));
        end
        yticklab{ncbarlabels} = sprintf('%s%1.2f',maxcbarlabel, ytick_level(ncbarlabels));
      else
        yticklab{1} = sprintf('%s%1.2e',mincbarlabel,ytick_level(1));
        for c = 2:ncbarlabels-1
          yticklab{c} = sprintf('%1.2e',ytick_level(c));
        end
        yticklab{ncbarlabels} = sprintf('%s%1.2e',maxcbarlabel, ytick_level(ncbarlabels));
      end
      
      if ~h % Generate a new figure
        h = figure('Name',fstr,'color','w','position',fpos);
        handles = [handles,h];
      else
        %set(0,'currentfigure',h);
      end
      
      if over && f ~= elem_ana_range(1)
        eltext = sprintf('%s, %s', eltext,fldnm{f});
      else
        eltext = fldnm{f};
        clear title_text
        %type = '';
      end
      
      if over
        colr = over_col;
        colr(:,over_logical(chanord(n),:) == 0) = 0; %mod(f-1,size(over_logical,1)),:)==0) = 0;
        %d(d <= r1) = r1; d(d >= r2) = r2; %reshape(thresh_d,size(d));
        dimg = repmat(mat2gray(d,[r1,r2]),[1,1,3]);
        nonzeroch = find(over_logical(chanord(n),:) == 1);
        colrsused(n,:) = over_logical(chanord(n),:);
        if n == 1
          dplot = zeros(size(d,1),size(d,2),3);
          PropertyNames{length(PropertyNames)+1} = 'h';
          PropertyVal{length(PropertyNames)+1} = h;
        end
        
        if length(elem_ana_range) == 1;
          dplot = d;
          if strmatch('Linear',ptag{pt})
            imagesc(x, y, dplot,'alphadata',trans_alpha);
            caxis(h,[r1,r2]);
            colormap(h,colr);
            %cbarh = colorbar('southoutside','fontsize',fntsz,'xtick',ytick_level,'xticklabel',yticklab);
            if cb
              switch cbl
                case 'southoutside'
                  cbarh = colorbar(cbl,'fontsize',fntsz,'xtick',xtick_level,'xticklabel',xticklab);
                case 'eastoutside'
                  cbarh = colorbar(cbl,'fontsize',fntsz,'ytick',ytick_level,'yticklabel',yticklab);
              end
            end
          end
          
          if strmatch('Log',ptag{pt})
            imagesc(x, y, log10(dplot),'alphadata',trans_alpha);
            caxis(h,log10([r1,r2]));
            colormap(h,colr);
            if cb; colorbar(cbl); end
          end
          
        else
          dplot(:,:,nonzeroch) = (dplot(:,:,nonzeroch) + dimg(:,:,nonzeroch));%./length(elem_ana_range);
          errori = find(dplot > 1);
          [dim1,dim2,dim3] = ind2sub(size(dplot),errori);
          d3 = 1;
          % This method of combining data keeps high pixels in each
          % individual color high, but reduces the magnitude of pixels
          % which are high in BOTH colors being shown. For example, a pixel
          % which is 100 for one element and 010 for another will not be
          % affected, but is a pixel is 110 for one element and 010 for
          % another, it will be reduced to 0.5,1,0. one problem with this
          % method is a color that is 0.5,0,0 for one element and 0,1,0
          % will have the same coloration
          while ~isempty(errori) && d3 <= 3
            evali = find(dim3 == d3);
            evalrows = dim1(evali);
            evalcols = dim2(evali);
            for e = 1:length(evalrows)
              dplot(evalrows(e),evalcols(e),:) = dplot(evalrows(e),evalcols(e),:)/dplot(evalrows(e),evalcols(e),d3);
            end
            errori = find(dplot > 1);
            [dim1,dim2,dim3] = ind2sub(size(dplot),errori);
            d3 = d3 + 1;
          end
          % This method scales all pixels down by a factor equal to the
          % maximum pixel value in any channel (RGB) in the entire dataset.
          % This method is not as great because it dim ALL pixel value for
          % both all channels
          %scale_factor = max(dplot(:));
          %dplot = dplot./scale_factor;
          if strmatch('Linear',ptag{pt})
            imagesc(x, y, dplot);
          end
          
          if strmatch('Log',ptag{pt})
            imagesc(x, y, log10(dplot));
          end
          
          if n == length(elem_ana_range);
            disptxt = sprintf(' %s',strrep(eltext,'_',''));
            commas = strfind(disptxt,',');
            commas = [1,commas,length(disptxt)+1];
            if ~ nomarks || legon
              for m = 1:n
                hold on;
                plot(x(1),y(1),'s','markerfacecolor',colrsused(m,:),...
                  'markeredgecolor',colrsused(m,:),'displayname',disptxt((commas(m)+1):(commas(m+1)-1)));
              end
            end
            legend('show')
            showtitle = 'h';
          end
        end
        
        if sc; axis image; end%set(gca,'dataaspectratio',[1,1,1]); end
        
      else
        if chanord == 8
          colr = jet(color_res);
        elseif chanord == 9
          colr = custom_colormaps('fire');
        elseif chanord == 10
          colr = custom_colormaps('gray');
        end
        dplot = d;
        if strmatch('Linear',ptag{pt})
          imagesc(x, y, dplot,'alphadata',trans_alpha);
          caxis([r1,r2]);
          colormap(colr);
          %cbarh = colorbar('southoutside','fontsize',fntsz,'xtick',ytick_level,'xticklabel',yticklab);
          if cb
            cbarh = colorbar(cbl,'fontsize',fntsz);
            %             switch cbl
            %               case 'southoutside'
            %                 cbarh = colorbar(cbl,'fontsize',fntsz,'xtick',ytick_level,'xticklabel',yticklab);
            %               case 'eastoutside'
            %                 cbarh = colorbar(cbl,'fontsize',fntsz,'ytick',ytick_level,'yticklabel',yticklab);
            %             end
          end
        end
        
        if strmatch('Log',ptag{pt})
          imagesc(x, y, log10(dplot),'alphadata',trans_alpha);
          caxis(log10([r1,r2]));
          colormap(colr);
          if cb; colorbar; end
        end
        
        if sc; set(gca,'dataaspectratio',[1,1,1]); end
      end
      xlabh = xlabel(xstr,'fontsize',fntsz); set(gca,'fontsize',fntsz);
      ylabh = ylabel(ystr,'fontsize',fntsz);
      %{
      ext_x = get(xlabh,'position');
      xLIM = xlim;
      txth = text('string',label(end-4:end),'position',[xLIM(2)*1.05,ext_x(2)],...
        'fontsize',fntsz)
      %}
      xLIM = xlim;
      yLIM = ylim;
      if over && length(elem_ana_range) == 1 && cb || ~over
        switch cbl
          case 'southoutside'
            txth = text('string',type(strfind(type,'['):strfind(type,']')),...
              'position',[xLIM(end),yLIM(end)],...
              'fontsize',fntsz,'VerticalAlignment','top','horizontalalignment','right');
          case 'westoutside'
            txth = text('string',type(strfind(type,'['):strfind(type,']')),...
              'position',[xLIM(2),mean(ylim)],...
              'fontsize',fntsz,'rotation',90,'VerticalAlignment','top');
        end
      end
      
      line_2text = strcat(thresh_str,line_str,sample_str,elim_str);
      if ~exist('title_text') || over
        if ~isempty(line_2text)
          title_text = {sprintf('%s %s %s %s %s Heatmap',data_str,type,ptag{pt},eltext,dmod),...
            line_2text};
          if ~isempty(strfind(title_text,'_'))
            underind = strfind(title_text,'_');
            for k = 1:length(title_text);
              tt_vec = 1:length(title_text{k});
              tt_vec = tt_vec(find(ismember(tt_vec,underind{k}) == 0));
              title_text{k} = title_text{k}(tt_vec);
            end
          end
        else
          title_text = sprintf('%s %s %s %s %s Heatmap',data_str,type,ptag{pt},eltext,dmod);
          if ~isempty(strfind(title_text,'_'))
            underind = strfind(title_text,'_');
            tt_vec = 1:length(title_text);
            tt_vec = tt_vec(find(ismember(tt_vec,underind) == 0));
            title_text = title_text(tt_vec);
          end
        end
      end
      
      if showtitle==1
        title(title_text);
      elseif strmatch(showtitle,'d')
        title(title_text(2:end))
      elseif strmatch(showtitle,'h')
        if ~iscell(title_text)
          title(title_text)
        else
          title(title_text(1))
        end
      end;
      %title({[sprintf('%s %s Linear Intensity [cps] Heatmap',data_str,fldnm{f})],...
      % [strcat(thresh_str,line_str,sample_str,elim_str)]},'fontsize',fntsz);
      
      if ~over
        validfldnm = regexprep(fldnm{f},'\s','_');
        dout = setfield(dout,validfldnm,d);
      else
        dout = dplot;
      end
    else
      if ~over
        validfldnm = regexprep(fldnm{f},'\s','_');
        dout = setfield(dout,validfldnm,d);
      end
    end
    
    if ~isempty(relhtxt)
      hold on;
      [line_types,~,linesind] = unique(relhtxt,'stable');
      nistind = find(cellfun(@(x) ~isempty(x),strfind(lower(line_types),'nist')));
      labelchangei = find(diff(linesind) ~= 0) + 1;
      starti = [1;labelchangei];
      endi = [labelchangei-1;length(relhtxt)];
      
      fntclr = 'w';
      bgclr = 'k';
      colrs = lines(length(starti));
      for i = 1:length(starti)
        if rem(i,2) == 0^(rem(nistind,2))
          %text(1.5,starti(i) + (endi(i) - starti(i) + 1)/2,strrep(relhtxt(starti(i)),'_','\_'),...
          text(N,starti(i),strrep(relhtxt(starti(i)),'_','\_'),...
            'fontsize',14,'HorizontalAlignment','right','BackgroundColor',bgclr,...
            'VerticalAlignment','top','Color',fntclr,'FontName','ariel');
          rectangle('position',[0.5,starti(i)-0.5,N+1,endi(i) - starti(i) + 1],...
            'edgecolor',colrs(i,:),'linewidth',2,'linestyle',':')
        end
      end
      svstr = 'labeled';
    end
    
    if showaqarrow
      hold on;
      lnwdth = 2;
      xdist = (xLIM(2) - xLIM(1))/50;
      ydist = (yLIM(2) - yLIM(1))/50;
      switch rotTag
        case -2
          x1 = xLIM(1) + xdist;
          y1 = yLIM(2) - ydist;
          x2 = xLIM(1) + xdist*4;
          y2 = y1;
          yt = y1 - ydist*2;
          xt = x1;
          %quiver(x1,y1,x2-x1,0,'color','r','linewidth',3)
          thal = 'left';
          tval = 'bottom';
          lhal = 'left';
          lval = 'bottom';
          lstr = '>';
          tstr = '^';
        case -1
          x1 = xLIM(2) - xdist;
          y1 = yLIM(1) + ydist;
          x2 = x1;
          y2 = yLIM(1) + ydist*4;
          xt = x1 - ydist*2;
          yt = y1;
          thal = 'right';
          tval = 'top';
          lstr = 'v';
          tstr = '<';
        case 0
          x1 = xLIM(1) + xdist;
          y1 = yLIM(1) + ydist;
          x2 = xLIM(1) + xdist*4;
          y2 = y1;
          xt = x1;
          yt = y1 + xdist*2;
          thal = 'left';
          tval = 'top';
          lstr = '>';
          tstr = 'v';
        case 1
          x1 = xLIM(1) + xdist;
          y1 = yLIM(2) - ydist;
          x2 = x1;
          y2 = y1 - ydist*4;
          xt = x1 + ydist*2;
          yt = y1;
          thal = 'left';
          tval = 'bottom';
          lstr = '^';
          tstr = '>';
        case 2
          x1 = xLIM(2) - xdist;
          y1 = yLIM(2) - ydist;
          x2 = x1 - xdist*4;
          y2 = y1;
          xt = x1;
          yt = y1 - ydist*2;
          thal = 'right';
          tval = 'bottom';
          lstr = '<';
          tstr = '^';
      end
      plot(x1,y1,'or','LineWidth',lnwdth,'HandleVisibility','off')
      plot(x2,y2,strcat(lstr,'r'),'LineWidth',lnwdth,'HandleVisibility','off')
      plot([x1,x2],[y1,y2],'r','LineWidth',lnwdth,'HandleVisibility','off')
      plot(xt,yt,strcat(tstr,'g'),'LineWidth',lnwdth/2,'HandleVisibility','off')
      plot([x1,xt],[y1,yt],'g','LineWidth',lnwdth/2,'HandleVisibility','off')
      text(xt,yt,'t','HorizontalAlignment',thal,'VerticalAlignment',tval,'color','g')
    end
    clear d
    if nomarks
      if ~legon
        legend off;
      end
      set(gca,'box','off','xtick',[],'ytick',[],'xcolor','w','ycolor','w',...
        'color','none','position',[0,0,1,1]);
      title('');
      if cb
        cbh = colorbar;
        cph = get(cbh,'position');
        set(cbh,'position',[min([cph(1),0.91]),0.1,0.02,0.8]);
      end
    end
  end
end

if sv && ~isempty(strmatch('Log',ptag{pt}))
  if over && sc
    svfstr = sprintf('spatial combo log %s %s',eltext,chanordstr);
  elseif over
    svfstr = sprintf('combo log %s %s',eltext,chanordstr);
  elseif sc
    svfstr = sprintf('spatial log');
  else
    svfstr = 'log';
  end
  
  if ~isempty(svstr)
    svfstr = sprintf('%s %s %s',svfstr,svstr,...
      strrep(type((strfind(type,'[')+1):(strfind(type,']')-1)),'/','per'));
  else
    svfstr = sprintf('%s %s',svfstr,...
      strrep(type((strfind(type,'[')+1):(strfind(type,']')-1)),'/','per'));
  end
  
  if dock
    dockf on all
  end
  save_open_figures(data,[],[],svfstr);
  if length(ptag) > 1; close all; end
elseif sv
  % ==== v8 and prior naming scheme
  if over && sc
    svfstr = sprintf('_%s_%s',eltext,chanordstr);%svfstr = sprintf('spatial combo %s %s',eltext,chanordstr);
  elseif over
    svfstr = sprintf('_%s_%s',eltext,chanordstr);%svfstr = sprintf('combo %s %s',eltext,chanordstr);
  elseif sc
    svfstr = '';%sprintf('spatial');
  else
    svfstr = '';
  end
  
  if ~isempty(svstr)
    svfstr = sprintf('%s %s',svstr,svfstr);%,...
    %strrep(type((strfind(type,'[')+1):(strfind(type,']')-1)),'/','per'));
  elseif ~isempty(svfstr)
    svfstr = sprintf('%s',svfstr);%,...
    %strrep(type((strfind(type,'[')+1):(strfind(type,']')-1)),'/','per'));
  else
    svfstr = '';%sprintf('%s',...
    %strrep(type((strfind(type,'[')+1):(strfind(type,']')-1)),'/','per'));
  end
  
  if dock
    dockf on all
  end
  save_open_figures(data,[],[],strrep(svfstr,' ',''),varargin{:});
end

varargout{1}.figh = handles;

if ~over && isstruct(data)
  if ~isfield(dout,'Time')
    dout = setfield(dout,'Time',getfield(data,'Time'));
  end
  
  if ~isfield(dout,'line_datevec')
    dout = setfield(dout,'line_datevec',getfield( data,'line_datevec'));
  end
  
  if ~isfield(dout,'dataset')
    dout = setfield(dout,'dataset',getfield( data,'dataset'));
  end
end
varargout{2} = dout;
varargout{3} = x;
varargout{4} = y;
if dock
  dockf on all
end