function [datafxed,fxedinfo] = fix_sngl_pxl(data,elem_ana_range,varargin);

%% Script:  [datafxed,fxedinfo] = fix_sngl_pxl(data,elem_ana_range,varargin);
% Description:  Fixes single pixel anomalies in data. First, looks at
% points withing an N-point window (where N is odd), and fits the data to a
% line using the first and last points in the window.  Any points above the
% line are flagged and their derivative of that point and the neighboring
% point is checked.  If the derivative is largely positive and then largely
% negative (with values within a point to point similarity, determined by
% the user), then that point is redefined as the mean of the neighboring
% points.
% TESTED for npnts = 1 and 2 on 041211\Abrain data with satisfactory output
% Example:
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotope analyzed
% varargin - 'PropertyName','PropertyValue'
%   • 'winlen':  numeric specifying window length [DEFAULT = 5]
%   • 'p2psim':  numeric fraction indicating the degree of point to point
%   similarity requires, used to find spikes in data [DEFAULT = 0.7]
%   • 'lines': row matrix indicating specific lines user would like to
%   analyze [DEFAULT = 'all']
%   • 'plot_lines': binary indicating whether to show plots of data in line with
%   informating indicating changed points and value of new points (1) or
%   not (0) [DEFAULT = 0]
%   • 'hist': binary indicating whether to show histogram of points number
%   of points changed for a given magnitude of change [DEFAULT = 0]
%   • 'hdim': numeric indicating whether to plot histograms in 2-D (one
%   histogram plot/isotope) or 3D (one histogram plot for all isotopes)
%   [DEFUALT = 3]
%   • 'npnts': number of consecutive high points to cut out [DEFAULT = 2]
%   • 'heatmap': binary indicating whether to show plots (1) or not (0) [DEFAULT
%   = 0]
%   • 'error_heatmap': binary indicating whether to plot the error heatmap
%   (1) or not (0) [DEFAULT = 0]
%   • 'save_data': binary indicating whether to save restructured data (1)
%   or not (0). Data are saved in the dataset directory with element names
%   (they overwrite existing raw element data csv files with the same name)
%   [DEFAULT = 0]
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  4 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  5 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  13 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     3 -
%  Attempting to correct redefinition of pixels to higher values
%  20 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     4


datafxed = [];
fxedinfo = [];
%fxed_all = [];
fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('winlen',PropertyNames)
  winlen = PropertyVal{strmatch('winlen',PropertyNames)};
else
  winlen = 5;
end

if strmatch('p2psim',PropertyNames) % point to point similarity. Can differ by this multiple.
  p2psim = PropertyVal{strmatch('p2psim',PropertyNames)};
else
  p2psim = 0.7;
end

if strmatch('lin',PropertyNames)
  lines = PropertyVal{strmatch('lin',PropertyNames)};
else
  lines = 1:size(getfield(data,fldnm{1}),1);
end

if strmatch('plot_l',PropertyNames)
  p = PropertyVal{strmatch('plot_l',PropertyNames)};
else
  p = 0;
end

if strmatch('error',PropertyNames)
  erhtm = PropertyVal{strmatch('error',PropertyNames)};
else
  erhtm = 0;
end

if strmatch('heat',PropertyNames)
  htm = PropertyVal{strmatch('heat',PropertyNames)};
else
  htm = 0;
end

if strmatch('hist',PropertyNames)
  h = PropertyVal{strmatch('hist',PropertyNames)};
  if strmatch('hdim',PropertyNames)
    hdim = PropertyVal{strmatch('hdim',PropertyNames)};;
  else
    hdim = 2;
  end
else
  h = 0;
end

if strmatch('np',PropertyNames)
  npnts = PropertyVal{strmatch('np',PropertyNames)};
else
  npnts = 2;
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 0;
end

for f = elem_ana_range
  info_struct = [];
  %out_data = [];
  %disp(fldnm{f});
  if isstruct(data)
    d = getfield(data,fldnm{f});
  else
    d = data;
  end
  row_col_all = [];
  row_col_10 = [];
  row_col_110 = [];
  row_col_100 = [];
  if isempty(strmatch(fldnm{f},{'Time','line_datevec','dataset'})) && f <= length(fldnm)
    row_col_ind = [];
    for l = lines
      lndat = d(l,:);
      lndiff = diff(lndat);
      pad_ln = [lndat(1).*ones(1,(winlen-1)/2),lndat,lndat(end).*ones(1,(winlen-1)/2)];
      pad_lndiff = [lndiff(1).*ones(1,(winlen+1)/2),lndiff,lndiff(end).*ones(1,(winlen-1)/2)];
      
      if p
        figure('color','w');
        plot([pad_ln;pad_lndiff]','linestyle','-','marker','.','linewidth',2);
        hold on;
        legend({'Line Intensity','Line Derivative'},'fontsize',fntsz);
        title(sprintf('%s 041211\\Abrain Line %d',fldnm{f},l),'fontsize',fntsz)
        set(gca,'fontsize',14)
        axis tight
      end
      
      ckind = [];
      fxind = [];
      fx_110 = [];
      fx_10 = [];
      fx_100 = [];
      %disp(l)
      for i = npnts:(length(pad_ln)-winlen-1)
        indrng = i:i+winlen-1;
        windat = pad_ln(indrng);
        m = (windat(end)-windat(1))/(winlen-1);
        b = (winlen*windat(1)-windat(end))/(winlen-1);
        lindat = m.*[1:winlen]+b;
        gtpnts = find(windat > lindat);
        gtind = gtpnts + i -1;
        for g = gtind(:)'
          if ~ismember(g,ckind)
            ckind = [ckind,g];
            der = pad_lndiff(g-npnts+1:g+1);%g-1:g+npnts-1);
            %% CASE 1 - Checks for a large spike in instensity followed by
            % a sharp drop. Sometimes, spikes will occur over npnts. Checks
            % if there is an npnt spike and then a sharp drop. Searches for
            % three point sequence of HIGH-HIGH-LOW
            if sum(der(1:npnts) > 0) == npnts && der(end) < 0
              if min(der(end-1),abs(der(end))) >= max(der(end-1),abs(der(end)))*p2psim
                fxind = [fxind,g];
                fx_10 = [fx_10,g];
                if p; plot(g,pad_ln(g),'*r'); end
                pad_ln(g) = (pad_ln(g-1)+pad_ln(g+1))/2;
                if p; plot(g,pad_ln(g),'*g'); end
              elseif min(sum(der(1:npnts))*0.8,abs(der(end))) >= max(sum(der(1:npnts))*0.8,abs(der(end)))*p2psim ...
                  && npnts > 1
                [pad_ln,fxind,fx_110] = fix_pnt(pad_ln,g,npnts,fxind,fx_110,p);
              end
              %% CASE 2 - checks if there's a single large spike and then a
              % large drop in intensity (single point). Searches for two
              % point sequence HIGH-LOW
            elseif der(end-1) > 0 && der(end) < 0
              if min(der(end-1),abs(der(end))) >= max(der(end-1),abs(der(end)))*p2psim
                fxind = [fxind,g];
                fx_10 = [fx_10,g];
                if p; plot(g,pad_ln(g),'*r'); end
                pad_ln(g) = (pad_ln(g-1)+pad_ln(g+1))/2;
                if p; plot(g,pad_ln(g),'*g'); end
              end
              %% CASE 3 - Checks if there's a large spike followed by two
              % neighboring sharp drops. Searching for three point sequence
              % HIGH-LOW-LOW
            elseif der(1) > 0 && sum(der(2:end) < 0) == length(der)-1
              if min(der(1),abs(sum(der(2:end))*0.8)) >= max(der(1),abs(sum(der(2:end))*0.8))*p2psim
                [pad_ln,fxind,fx_100] = fix_pnt(pad_ln,g,npnts,fxind,fx_100,p);
              end
            end
          end
        end
      end
      d_corrected(l,:) = pad_ln(1+floor(winlen/2):end-floor(winlen/2));
      row_col_all = [row_col_all,[l.*(ones(1,length(fxind)));fxind-floor(winlen/2)]];
      row_col_10 = [row_col_10,[l.*(ones(1,length(fx_10)));fx_10-floor(winlen/2)]];
      row_col_110 = [row_col_110,[l.*(ones(1,length(fx_110)));fx_110-floor(winlen/2)]];
      row_col_100 = [row_col_100,[l.*(ones(1,length(fx_100)));fx_100-floor(winlen/2)]];
      if p; plot(pad_ln,'m'); end;
    end
    
    
    %% Quantify amount of denoising performed in the image
    dfrac = (d-d_corrected)./d;
    if ~isempty(row_col_all)
      dfracvec = dfrac((row_col_all(2,:)-1).*size(d,1)+row_col_all(1,:));
      d10 = dfrac((row_col_10(2,:)-1).*size(d,1)+row_col_10(1,:));
      d100 = dfrac((row_col_100(2,:)-1).*size(d,1)+row_col_100(1,:));
      d110 = dfrac((row_col_110(2,:)-1).*size(d,1)+row_col_110(1,:));
    else
      dfracvec = [];
      d10 = [];
      d100 = [];
      d110 = [];
    end
    
    %% Shape output variables
    info_struct = setfield(info_struct,'fxed_all',row_col_all);
    info_struct = setfield(info_struct,'fxed_10',row_col_10);
    info_struct = setfield(info_struct,'fxed_110',row_col_110);
    info_struct = setfield(info_struct,'fxed_100',row_col_100);
    info_struct = setfield(info_struct,'dfrac',dfrac);
    info_struct = setfield(info_struct,'dfracsum',sum(abs(dfracvec)));
    info_struct = setfield(info_struct,'perc_pnts',(size(row_col_all,2)/(size(d,1)*size(d,2))).*100);
    info_struct = setfield(info_struct,'perc_10',(length(d10)/length(dfracvec))*100);
    info_struct = setfield(info_struct,'perc_100',(length(d100)/length(dfracvec))*100);
    info_struct = setfield(info_struct,'perc_110',(length(d110)/length(dfracvec))*100);
    info_struct = setfield(info_struct,'frac_mean', mean(dfracvec));
    info_struct = setfield(info_struct,'frac_npos', length(dfracvec(dfracvec>0)));
    info_struct = setfield(info_struct,'frac_nneg', length(dfracvec(dfracvec<0)));
    info_struct = setfield(info_struct,'frac_nzero', length(dfracvec(dfracvec==0)));
    datafxed = setfield(datafxed,fldnm{f},d_corrected);
    fxedinfo = setfield(fxedinfo,fldnm{f},info_struct);
    %info_struct = setfield(info_struct,'data', length(dfracvec(dfracvec==0)));
    %fxed_all = setfield(fxed_all,fldnm{f},info_struct);
    if htm; plot_heatmap(d_corrected,sprintf('%s denoised',fldnm{f}),'plot','lin','thresh','auto',varargin{:}); end
    if erhtm; plot_heatmap(dfrac,sprintf('%s error',fldnm{f}),'plot','lin','thresh','auto',varargin{:}); end
  elseif strmatch(fldnm{f},'Time')
    datafxed = setfield(datafxed,'Time',data.Time);
    fxedinfo = setfield(fxedinfo,'Time',data.Time);
    %out_data = setfield(out_data,'data',data.Time);
    %fxed_all = setfield(fxed_all,'Time',out_data);
  elseif strmatch(fldnm{f},'line_datevec')
    datafxed = setfield(datafxed,'line_datevec',data.line_datevec);
    fxedinfo = setfield(fxedinfo,'line_datevec',data.line_datevec);
  elseif strmatch(fldnm{f},'dataset')
    datafxed = setfield(datafxed,'dataset',data.dataset);
    fxedinfo = setfield(fxedinfo,'dataset',data.dataset);
  end
end
    
datafxed = setfield(datafxed,'Time',data.Time);
fxedinfo = setfield(fxedinfo,'Time',data.Time);

datafxed = setfield(datafxed,'line_datevec',data.line_datevec);
fxedinfo = setfield(fxedinfo,'line_datevec',data.line_datevec);

datafxed = setfield(datafxed,'dataset',data.dataset);
fxedinfo = setfield(fxedinfo,'dataset',data.dataset);

% if htm; plot_heatmap(datafxed,sprintf('%s denoised',fldnm{f}),'plot','lin','thresh','auto'); end
% if erhtm; plot_heatmap(dfrac,sprintf('%s error',fldnm{f}),'plot','lin','thresh','auto'); end
if h
  plot_fix_sngl_pxl(fxedinfo,elem_ana_range,varargin{:});
end

if svdat
  %loaddir; home_dir = strcat(ladatapath,datafxed.dataset,'\');
  currentfldr = cd;
  Dpos = strfind(currentfldr,'Documents');
  
  home_dir = strcat(currentfldr(1:Dpos-1),'Documents\MADLab Data\Laser Ablation\',datafxed.dataset,'\');
  allfiles = ls;
  files = cellstr(allfiles);
  
  for i = 1:length(fldnm)
    if isempty(strmatch(fldnm{i},{'dataset'}))
      csvname = strcat(home_dir,fldnm{i},'.csv');
      csvwrite(csvname,getfield(datafxed,fldnm{i}));
    end
  end
end

%% SUBFUNCTIONS-----------------------------------------------
  function [pad_ln,fxind,fxind_spec] = fix_pnt(pad_ln,g,npnts,fxind,fxind_spec,p);
    if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*r'); end
    y1 = pad_ln(g-npnts); y2 = pad_ln(g+1); x1 = g-npnts; x2 = g+1;
    m = (y2-y1)/(x2-x1);
    b = (x2*y1-x1*y2)/(x2-x1);
    new_pnts = m.*[g-1,g]+b;
    nznew_pnt = find((pad_ln([g-1,g])-new_pnts) > 0);
    i = [g-1,g];
    pad_ln(i(nznew_pnt)) = new_pnts(nznew_pnt);
    fxind = [fxind,i(nznew_pnt)];
    fxind_spec = [fxind_spec,i(nznew_pnt)];
    if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*g'); end
  end
end

