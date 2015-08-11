function varargout = plot_ivisimg(data,img_ana_ind,varargin)

%% Script: [varargout{figure handle, plotted data, x coords, y coords}]=
%   plot_heatmap(data,elem_ana_range,varargin)
% Description:
% Example:
% Required Function: interp_data_elemrng,
%
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% img_ana_ind: 1 x <# photos> vector which corresponds to the field
%       number in the data structure. This variable can also be a cell of
%       strings where the string specifies the images the user would like
%       to plot
%       1 = 'fluorescentreference' ('f')
%       2 = 'luminescent' ('l')
%       3 = 'photograph' ('p')
%       4 = 'readbiasonly' ('r')
% varargin - 'PropertyName','PropertyValue'
%   • 'color_res':  numeric specifying color resolution of the heatmap
%   (maximum 256) [DEFAULT = 256]
%   • 'thresh_data':  1 x 2 matrix with minimum threshold value as first
%   element and maximum threshold value as second. Can also use string
%   'auto' to cut off bottom 0.02% and top 0.5% of data points
%   • 'plot_type': string - 'log','linear','both' [DEFUALT: 'lin']
%   • 'scale': 1 or 0 - uses spatial scaling rather than line # and sample #
%   (1) or not (0) [DEFAULT = 0]
%   • 'overlay': binary indicating whether to overlay input isotopes (1) or
%   not (0) [DEFAULT = 0]
%   • 'handle': figure handle user would like to generate the images on. (0)
%   to generate a new figure. [DEFAULT = 0]
%   • 'type': string indicating whether data is given in insensity (counts
%   per second) -- 'i'-- or absolute concentration (parts per million)
%   --'c'-- [DEFAULT = 'i']
%   • 'imgnum': vector indicating the image number that should be plotted
%   [DEFAULT = all images in data.imfield.image]
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
%   • 'dockall': binary indicationg whether to dock all figures (1) or not
%   (0) [DEFAULT = 0]
%   • 'figpos': 4 x 1 vector indicating the figure position [DEFAULT = []]
%   • 'showcbar': binary indicating whether to show the colorbar (1) or
%   not (0) [DEFAULT = 1]
%
% OUTPUTS ---------------------------------------------------------------
% varargout
%   {1}: Figure handle
%
%  Date           Author            E-mail                      Version
%  26 Nov  2014   Amanda Balderrama amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

fldnm = {'fluorescentreference','luminescent','photograph','readbiasonly'};
if ~isempty(img_ana_ind) && isstruct(data)
  if isnumeric(img_ana_ind)
    fldnm = fldnm(img_ana_ind);
  elseif isstr(img_ana_ind) || iscell(img_ana_ind)
    if iscell(img_ana_ind)
      img_ana_ind = cellfun(@(x) x(1),img_ana_ind);
    end
    
    img_ana_str = img_ana_ind;
    img_ana_ind = zeros(1,length(img_ana_str));
    for i = 1:length(img_ana_ind)
      iout = strfind(fldnm,img_ana_str(i));
      inonemp = nonemptystrind(fldnm,img_ana_str(i));
      firsti10 = cellfun(@(x) any(x == 1),iout(inonemp));
      img_ana_ind(i) = inonemp(firsti10);
    end
  end
elseif isstruct(data)
  img_ana_ind = 1:length(fldnm);
else
  img_ana_ind = 1;
  fldnm = 'din';
end


%% Initialize variables that deal with scaling the image (converting pixels --> cm)
if strmatch('scale',PropertyNames)
  sc = PropertyVal{strmatch('scale',PropertyNames)};;
else
  sc = 0;
end

if strmatch('showtitle',PropertyNames)
  showtitle = PropertyVal{strmatch('showtitle',PropertyNames)};;
else
  showtitle = 1;
end

if strmatch('title_txt',PropertyNames)
  title_text = PropertyVal{strmatch('title_txt',PropertyNames)};;
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
  fpos = [360,502,560,420];
end

if strmatch('showcbar',PropertyNames)
  cb = PropertyVal{strmatch('showcbar',PropertyNames)};
else
  cb = 1;
end

%% Color parameters
if strmatch('color',PropertyNames)
  color_res = PropertyVal{strmatch('color',PropertyNames)};
  if color_res > 256
    color_res = 256;
  end
else
  color_res = 256;
end

%% Overlay parameters
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

% if strmatch('type',PropertyNames)
%   type = PropertyVal{strmatch('type',PropertyNames)};
% else
%   type = 'i'; %intensity [cps] versus concentration [ppm]
% end
%
% if strmatch('c',type)
%   type = 'Concentration [ppm]';
% else
%   type = 'Intensity [cps]';
% end

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
  chanordstr = chanordstr(1:length(img_ana_ind));
  chanord = 1:7;
else
  chanord = 8;
end

n = 0;
for pt = 1:length(ptag)
  for f = img_ana_ind
    n = n + 1;
    if isstruct(data)
      imgstruct = getfield(data,fldnm{f});
      imfields = fieldnames(imgstruct);
      
      if isfield(imgstruct,'ImageUnits');
        IUcell = getfield(imgstruct,'ImageUnits');
        type = IUcell{1};
      else type = 'counts';
      end
      
      dcell = getfield(imgstruct,'image');
      switch f
        case 1 %fluorescentreference
          relflds = {'BinningFactor','Emissionfilter','Excitationfilter',...
            'FieldofView','Fluorescencelevel','Fluorescencelevelstring',...
            'ImageSubtype','ImageUnits','LuminescentExposureSeconds',...
            'Subjectsize','fNumber'};
        case 2 %luminescent
          relflds = {'BinningFactor','Emissionfilter','Excitationfilter',...
            'FieldofView','Fluorescencelevel','Fluorescencelevelstring',...
            'ImageSubtype','ImageUnits','LuminescentExposureSeconds',...
            'Subjectsize','fNumber'};
        case 3 %photograph
          relflds = {'BinningFactor','Emissionfilter','ExposureTimeSec',...
            'FieldofView',...
            'ImageUnits','Objectradius',...
            'Subjectsize','fNumber'};
        case 4 %readbiasonly
          relflds = {'BinningFactor','DataMultiplier','BackgroundExposureSeconds'};
      end
    else
      dcell = {d};
    end
    
    if strmatch('imgnum',PropertyNames)
      dimgv = PropertyVal{strmatch('imgnum',PropertyNames)};
    else
      dimgv = 1:length(dcell);
    end
    nimgs = length(dimgv);
    dimgsz = size(dcell{1});
    dc2m = double(cell2mat(dcell(dimgv)));
    
    thresh_str = '';
    if strmatch('thresh',PropertyNames)
      td = PropertyVal{strmatch('thresh',PropertyNames)};
      thresh_data = td(min(size(td,1),n),:);
      thresh_d = dc2m(:);
      if ~isempty(strmatch('auto',thresh_data)) || (sum(thresh_data - 10) < 0)...
          || (ischar(td) && isempty(td))
        thresh_data = auto_thresh(dc2m,[0.02,0.5]./nimgs,[]); % if 'auto', thresholds are at [0.02,0.5]
      else
        thresh_d(1,1) = thresh_data(1);
        thresh_d(2,1) = thresh_data(2);
      end
      ind_gtT = find(thresh_d > max(thresh_data));
      ind_ltT = find(thresh_d < min(thresh_data));
      thresh_d(ind_gtT) = max(thresh_data);
      thresh_d(ind_ltT) = min(thresh_data);
      dth = reshape(thresh_d,[dimgsz(1),dimgsz(2),nimgs]);
      thresh_str = sprintf('THRESHOLDS: lower = %10.3g; upper =  %10.3g; ',min(thresh_data),max(thresh_data));
      maxcbarlabel = '> ';
      if isempty(ind_ltT) || ~isempty(strfind(type,'ppm'))
        mincbarlabel = '';
      else
        mincbarlabel = '< ';
      end
    else
      maxcbarlabel = '';
      mincbarlabel = '';
      dth = reshape(dc2m,[dimgsz(1),dimgsz(2),nimgs])
    end
    
    if isfield(imgstruct,'FieldofView')
      fov = getfield(imgstruct,'FieldofView');
    else
      fov = repmat(dimgsz(1),1,nimgs);
    end
    
    if isfield(imgstruct,'BinningFactor')
      binningfact = getfield(imgstruct,'BinningFactor');
    else
      binningfact = ones(1,nimgs);
    end
    
    if isfield(imgstruct,'Excitationfilter')
      exfilt = getfield(imgstruct,'Excitationfilter');
    else
      exfilt = nan(1,nimgs);
    end
    
    if isfield(imgstruct,'Emissionfilter')
      emfilt = getfield(imgstruct,'Emissionfilter');
    else
      emfilt = nan(1,nimgs);
    end
    
    expose_field = nonemptystrind(imfields,'Exposure');
    second_field = nonemptystrind(lower(imfields),'sec');
    bothexsec = ismember(expose_field,second_field);
    if any(bothexsec)
      fldind = expose_field( bothexsec);
      exposevec = getfield(imgstruct,imfields{fldind});
    else
      exposevec = nan(1,nimgs);
    end
    
    if isfield(imgstruct,'fNumber')
      fnum = getfield(imgstruct,'fNumber');
    else
      fnum = nan(1,nimgs);
    end
    
    for dimg = 1:nimgs%dimgv(:)'
      d = dth(:,:,dimg);
      %% Figure handle
      if strmatch('h',PropertyNames)
        h = PropertyVal{strmatch('h',PropertyNames)};
      else
        h = 0;
      end
      
      if isempty(strmatch(fldnm{f},skip_fields))
        data_str = '';
        if isfield(data,'dataset')
          data_str = getfield(data,'dataset');
        end
        
        if iscell(exfilt) || ~isnan(exfilt(dimg))
          if iscell(exfilt)
            exstr = sprintf('Ex = %s,',exfilt{dimg});
          else
            exstr = sprintf('Ex = %d,',exfilt(dimg));
          end
        else exstr = '';
        end
        
        if iscell(emfilt) || ~isnan(emfilt(dimg))
          if iscell(emfilt)
            emstr = sprintf('Em = %s',emfilt{dimg});
          else
            emstr = sprintf('Em = %d',emfilt(dimg));
          end
        else emstr = '';
        end
        
        if ~isnan(exposevec(dimg))
          exposestr = sprintf('Exp(s) = %1.2f,',exposevec(dimg));
        else exposestr = '';
        end
        
        if ~isnan(fnum(dimg))
          fnumstr = sprintf('f# = %d,',fnum(dimg));
        else fnumstr = ''
        end
        
        if ~isnan(binningfact(dimg))
          binstr = sprintf('BinFact = %d,',binningfact(dimg));
        else binstr = ''
        end
        
        data_descrp = sprintf('%s %s %s %s %s', binstr, exposestr, fnumstr, exstr,emstr);
        
        if ~isempty(exstr) && ~iscell(exfilt)
          fstrmod = sprintf('_ex%d',exfilt(dimg));
        else
          fstrmod = '';
        end
        
        if over
          fstr = 'Heat'; % figure string
        else
          if isempty(dmod)
            fstr = sprintf('%s%d%s Heat', fldnm{f},dimgv(dimg),fstrmod); % figure string
          else
            fstr = sprintf('%s%d%s %s Heat',fldnm{f},dimgv(dimg),fstrmod,dmod); % figure string
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
        
        nonbinnedsz = size(d).*binningfact(dimg);
        if sc
          cmperpxl = fov(dimg)/nonbinnedsz(1);
          unscaled = 0:cmperpxl:fov(dimg);
          y = unscaled(1:binningfact(dimg):end);%[0,fov(dimg)];
          x = unscaled(1:binningfact(dimg):end);%[0,fov(dimg)];
          xstr = 'Distance [cm]';
          ystr = 'Distance [cm]';
        else
          x = 1:size(d,2);
          y = 1:size(d,1);
          xstr = '';
          ystr = '';
        end
        
        if over && f ~= img_ana_ind(1)
          eltext = sprintf('%s, %s', eltext,fldnm{f});
        else
          eltext = fldnm{f};
          clear title_text
          %type = '';
        end
        
        if over
          colr = over_col;
          colr(:,over_logical(chanord(n),:) == 0) = 0; %mod(f-1,size(over_logical,1)),:)==0) = 0;
          dimg = repmat(mat2gray(d),[1,1,3]);
          nonzeroch = find(over_logical(chanord(n),:) == 1);
          colrsused(n,:) = over_logical(chanord(n),:);
          if n == 1
            dplot = zeros(size(d,1),size(d,2),3);
            PropertyNames{length(PropertyNames)+1} = 'h';
            PropertyVal{length(PropertyNames)+1} = h;
          end
          
          if length(img_ana_ind) == 1;
            dplot = d;
            if strmatch('Linear',ptag{pt})
              imagesc(x, y, dplot,'cdata',[r1,r2],'alphadata',trans_alpha);
              colormap(h,colr);
              %cbarh = colorbar('southoutside','fontsize',fntsz,'xtick',ytick_level,'xticklabel',yticklab);
              if cb
                cbarh = colorbar('fontsize',fntsz,'ytick',ytick_level,'yticklabel',yticklab);
              end
            end
            
            if strmatch('Log',ptag{pt})
              imagesc(x, y, log10(dplot),'cdata',log10([r1,r2]),'alphadata',trans_alpha);
              colormap(colr);
              if cb; colorbar; end
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
            
            if n == length(img_ana_ind);
              disptxt = sprintf(' %s',eltext);
              commas = strfind(disptxt,',');
              commas = [1,commas,length(disptxt)+1];
              if ~ nomarks
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
            imagesc(x, y, dplot,'cdata',[r1,r2],'alphadata',trans_alpha);
            colormap(colr);
            %cbarh = colorbar('southoutside','fontsize',fntsz,'xtick',ytick_level,'xticklabel',yticklab);
            if cb
              cbarh = colorbar('fontsize',fntsz,'ytick',ytick_level,'yticklabel',yticklab);
            end
          end
          
          if strmatch('Log',ptag{pt})
            imagesc(x, y, log10(dplot),'cdata',log10([r1,r2]),'alphadata',trans_alpha);
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
        if over && length(img_ana_ind) == 1 && cb || ~over
          txth = text('string',type,'position',[xLIM(2),mean(ylim)],...
            'fontsize',fntsz,'rotation',90,'VerticalAlignment','top');
        end
        
        line_2text = strcat(thresh_str);
        if ~exist('title_text') || over
          if ~isempty(line_2text)
            title_text = {sprintf('%s %s %s %s %s Heatmap',data_str,type,ptag{pt},eltext,dmod),...
              data_descrp, line_2text};
          else
            title_text = {sprintf('%s %s %s %s %s Heatmap',data_str,type,ptag{pt},eltext,dmod),...
              data_descrp};
          end
          title_text = strrep(title_text,'_','\_');
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
          if ~isempty(exstr) && ~iscell(exfilt)
            validfldnm = sprintf('%s%d_ex%d',validfldnm,dimgv(dimg),exfilt(dimg));
          else
            validfldnm = sprintf('%s%d',validfldnm,dimgv(dimg));
          end
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
      clear d
    end
    
    if nomarks
      legend off;
      set(gca,'box','off','xtick',[],'ytick',[],'xcolor','w','ycolor','w',...
        'color','none','position',[0,0,1,1]);
      title('');
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
        svfstr = sprintf('%s %s %s',svfstr,svstr,type(end-3:end-1));
      else
        svfstr = sprintf('%s %s',svfstr,type(end-3:end-1));
      end
      
      if dock
        dockf on all
      end
      save_open_figures(data,[],[],svfstr);
      if length(ptag) > 1; close all; end
    elseif sv
      if over && sc
        svfstr = sprintf('spatial combo %s %s',eltext,chanordstr);
      elseif over
        svfstr = sprintf('combo %s %s',eltext,chanordstr);
      elseif sc
        svfstr = sprintf('spatial');
      else
        svfstr = '';
      end
      
      if ~isempty(svstr)
        svfstr = sprintf('%s %s %s',svfstr,svstr,type(end-3:end-1));
      elseif ~isempty(svfstr)
        svfstr = sprintf('%s %s',svfstr,type);
      else
        svfstr = sprintf('%s',type);
      end
      
      if dock
        dockf on all
      end
      save_open_figures(data,[],[],svfstr,varargin{:});
    end
  end
end

varargout{1}.figh = handles;

if ~over && isstruct(data)
  if ~isfield(dout,'origseq')
    dout = setfield(dout,'origseq',getfield( data,'origseq'));
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
end