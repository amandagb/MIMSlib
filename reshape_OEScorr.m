function [d,varargout] = reshape_OEScorr(din,varargin)
%% Script: d = reshape_OEScorr(dstand,varargin)
% [ds,reshape_info] = reshape_OEScorr(d,'lines',#,'temp_el','str', 'plot',0,'save_data',0,'save_heat',0,'top_perc',80,'bottom_perc',10);
% Description: Takes output from readfilesOES and restructures the data
% set. Note that results are quite sensitive to top and bottom percentages
% Example: [ds,reshape_info] = reshape_OEScorr(d,'lines',10,'temp_el','C_', 'plot',0,'save_data',0,'save_heat',0,'top_perc',80,'bottom_perc',10);
% Required Functions: find_linespersamp, find_rise
% Subfunctions: xcorr_align
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   • 'lines': integer indicating number of lines per sample [DEFAULT =
%   find_linespersamp function]
%   • 'temp_el': string indicate the template element to be used. The string
%   must uniquely identify the element and must have the first letter
%   capitalized. [DEFAULT = 'C_']
%   • 'plot': binary indicating whether to generate debugging plots (1) or
%   not (0) [DEFAULT = 0]
%   • 'save_heat': binary indicating whether to save heatmaps (1) or not (0)
%   [DEFAULT = 0]
%   • 'save_data': binary indicating whether to save restructured data (1)
%   or not (0). Data are saved in the dataset directory with element names
%   (they overwrite existing raw element data csv files with the same name)
%   [DEFAULT = 0]
%   • 'eval_samp_num': evaluates specific sample numbers [DEFAULT = evaluate
%   all sample numbers]
%   • 'temp_lines': determine which lines to use for the template
%   • 'keep_lines': vector containing lines that should be kept of final
%   image [DEFAULT = 'all']
%   • 'nalign': indicates the number of points that are allowed to vary
%   between the rises of subsequent lines [DEFAULT = 20]
%   • 'top_perc': integer indicating the percentage of data that should be
%   above the upper threshold. [DEFAULT = automated selection*--still in
%   debugging stages as of v5.4, 23 May 2012--use 20120320_OESCugrid40umps
%   for validation. The automated selection is done in lines 226-248]
%   • 'bottom_perc': integer indicating the percentage of data that should be
%   below the lower threshold. [DEFAULT = automated selection*--still in
%   debugging stages as of v5.4, 23 May 2012--use 20120320_OESCugrid40umps
%   for validation. The automated selection is done in lines 202-224]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  11 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  13 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Attempting to make template matching on all subsequent portions of
%     the data rather than fixed matching and template lines
%  14 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Attempts to align the derivative (either first or second) to a fixed
%     position so that left (rising) edge of the image is a straight line.
%  14 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.2
%     Trying to determine a better way for finding the initial rise of the
%     signal
%  20 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.3
%     Improving the rise point identification -- still several problems
%  21 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.4
%     Added rules for number of points before and after that must be less
%     than (if before) or greater than (if after) the current rise point
%  23 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.5
%     Adding functionality to view outputs for debugging as well as p
%  28 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.51
%     See line 388
%  28 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Will combine method used in reshape_OESstructure with line by line
%     correlation method -- 29 Mar: limited success, perhaps due to error
%     metric or length/uniformity of template from line to line
%  29 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5
%     See Mar 19 2012 Onenote notes about proposed solutions --
%     implementing ideal #2
%  31 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.2
%     Attempting to work with median filtered lines
%  10 Apr  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.21
%     Debugging script -- DEBUGGING PNTS:
%  09 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.3
%     Shortened the script so that processing occurs for each line
%     regardless of whether it's the beginning of a sequence or the end
%  14 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.4
%     Automating number of lines and top and bottom threshold
%     ** NEEDS TO BE DONE, good testing dataset(s) -- 20120305_OESAugrid,
%     Attempting to deal with the fact that the first few lines tend to
%     have significantly lower signal than subsequent lines (discovered
%     when debugging 20120501_OESpancreasR)
%  23 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.5
%     When working with 20120522_OESmbrain, it was found that while a
%     majority of the lines align correctly, only a few may be problematic.
%     It was found that misaligned lines can be identified by the parameter
%     "align_ind"'s 1st column == 0.
%  13 Dec  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.6
%  8  Jan  2014   Amanda Gaudreau   amanda.gaudreau@gmail.com     5.7
%     Working to make script functional for '20100312_kidney' dataset


close all

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

fldnm = fieldnames(din);

if strmatch('temp_el',PropertyNames)
  temp_el = PropertyVal{strmatch('temp_el',PropertyNames)};
else
  temp_el = 'C_';
end

if strmatch('lines',PropertyNames)
  linespersamp = PropertyVal{strmatch('lines',PropertyNames)};
else
  % Automated way to determine number of lines per sample
  data = getfield(din,fldnm{interp_data_elemrng(din,temp_el)});
  linespersamp = find_linespersamp(data);
end

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end

if strmatch('save_heat',PropertyNames)
  svheat = PropertyVal{strmatch('save_heat',PropertyNames)};
else
  svheat = 0;
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 0;
end

if strmatch('eval_sam',PropertyNames)
  eval_samp_num = PropertyVal{strmatch('eval_sam',PropertyNames)};
else
  eval_samp_num = 1:size(din.Time,1);
end

if strmatch('temp_lines',PropertyNames)
  temp_lines = PropertyVal{strmatch('temp_lines',PropertyNames)};
else
  temp_lines = 1:size(din.Time,1);
end

if strmatch('keep_lines',PropertyNames)
  keep_lines = PropertyVal{strmatch('keep_lines',PropertyNames)};
else
  keep_lines = 'a';
end

if strmatch('top_perc',PropertyNames)
  top_perc = PropertyVal{strmatch('top_perc',PropertyNames)};
else
  top_perc = [];
end

if strmatch('bottom_perc',PropertyNames)
  bottom_perc = PropertyVal{strmatch('bottom_perc',PropertyNames)};
else
  bottom_perc = [];
end

pntspersamp = size(din.Time,2) - 1;
if strmatch('pnts',PropertyNames)
  pntsperline = PropertyVal{strmatch('pnts',PropertyNames)};
else
  pntsperline = round(pntspersamp/linespersamp);
end

if strmatch('nalign',PropertyNames)
  nalign = PropertyVal{strmatch('nalign',PropertyNames)};
else
  nalign = 20;
end

if strmatch('xcorr',PropertyNames)
  xc = PropertyVal{strmatch('xcorr',PropertyNames)};
else
  xc = 0;
end

if strmatch('xcorrline',PropertyNames)
  linexc = PropertyVal{strmatch('xcorr',PropertyNames)};
else
  linexc = [];
end

d = [];
flag_lines = [];
himg = figure('position',[ 25   475   560   420]);%390   573   560   420]);%[403   246   560   420]);%
hlines = figure('position',[ 605   475   560   420]);%390    60   560   420]);%[403  -268   560   420]);%

if any(size(din.Time) == 1) || din.Time(1,1) == din.Time(2,1) % This condition indicates that the dataset has not previously been reshaped
  %% Generate Pattern using template element
  for f = interp_data_elemrng(din,temp_el)
    C = getfield(din,fldnm{f});
    if isempty(strmatch(fldnm{f},{'Time','line_datevec','dataset'}))
      C = C(:,1:pntspersamp);
      N = size(C,2);
      L = size(C,1);
      if any(size(din.Time) == 1)
        Cf = C;
      else
        Cf = medfilt1(C,5,L,2);
      end
      s = sort(Cf(:));
      Ncomp = linespersamp*L;
      
      % Determines the appropriate top and bottom thresholds that give the number
      % of lines per sample needed
      %{%
      eval_samp_num = find(sum(Cf > mean(C(:)),2));
      skipped_samples = find(ismember(1:L,eval_samp_num)==0);
      Cfvalid = Cf(eval_samp_num,:);
      if isempty(bottom_perc)
        ind_LT_mean = find(s <= mean(C(:)));
        bottom_perc0 = round(length(ind_LT_mean)/length(s)*100);
        bottom_perc = bottom_perc0;
        nLNobj = zeros(bottom_perc,length(eval_samp_num));
        while bottom_perc
          Nbot = round(length(C(:))*bottom_perc/100);  % number of points which give bottom percentage of data
          bot_thresh = s(Nbot+1);
          tfC = Cfvalid < bot_thresh;
          for j = 1:length(eval_samp_num)
            CC = bwconncomp(tfC(j,:));
            nLNobj(abs(-bottom_perc+bottom_perc0+1),j) = CC.NumObjects;
          end
          bottom_perc = bottom_perc - 1;
        end
        sum_closest = sum(abs(nLNobj-linespersamp-1),2);
        [v,bot_ind] = min(sum_closest);
        rep_min = (sum_closest == v);
        CC = bwconncomp(rep_min);
        obj_num = cellfun(@(x) ismember(bot_ind,x),CC.PixelIdxList);
        bot_ind = CC.PixelIdxList{find(obj_num==1)}(end);
        bottom_perc = bottom_perc0 + 1 - bot_ind;
      end
      
      if isempty(top_perc)
        ind_GT_mean = find(s > mean(C(:)));
        top_perc0 = round(length(ind_GT_mean)/length(s)*100);
        top_perc = top_perc0;
        nLNobj = zeros(100-top_perc-1,length(eval_samp_num));
        while (top_perc < 100-bottom_perc)
          Ntop = round(length(C(:))*top_perc/100);  % number of points which give bottom percentage of data
          top_thresh = s(length(s)-Ntop);
          tfC = Cfvalid > top_thresh;
          for j = 1:length(eval_samp_num)
            CC = bwconncomp(tfC(j,:));
            nLNobj(abs(top_perc-top_perc0+1),j) = CC.NumObjects;
          end
          top_perc = top_perc + 1;
        end
        sum_closest = sum(abs(nLNobj-linespersamp),2);
        [v,top_ind] = min(sum_closest);
        %rep_max = (sum_closest == v);
        %CC = bwconncomp(rep_max);
        %obj_num = cellfun(@(x) ismember(top_ind,x),CC.PixelIdxList);
        %top_ind = CC.PixelIdxList{obj_num}(end);
        top_perc = top_perc0 - 1 + top_ind;
      end
      %}
      disp(sprintf('Top percent = %d, Bottom percent: %d',top_perc,bottom_perc));
      Ntop = round(length(C(:))*top_perc/100);  % number of points which give top percentage of data
      Nbot = round(length(C(:))*bottom_perc/100);  % number of points which give bottom percentage of data
      top_thresh = s(length(s)-Ntop);
      bot_thresh = s(Nbot+1);
      
      ind = zeros(L*linespersamp,2);
      line_buffer = round(pntsperline/10); % used as a line buffer
      starti = 1;
      endi = pntsperline + line_buffer;
      pntsperline = endi;
      Y = zeros(L*linespersamp,pntsperline);
      Yxcorr = zeros(L*linespersamp,pntsperline);
      align_ind = zeros(L*linespersamp,2);
      rise_diffv = zeros(L*linespersamp,2);
      alskipped = []; %p = 0;
      yf = Cf(eval_samp_num(1),starti:endi);
      y = C(eval_samp_num(1),starti:endi);
      valid_lines = [];
      for l = eval_samp_num(:)'
        valid_lines = [valid_lines,(l-1)*linespersamp + 1:l*linespersamp];
        S = C(l,:);
        Sf = Cf(l,:);
        if sum(Sf >= top_thresh) > 0
          for k = 0:linespersamp-1
            line_num = (l-1)*linespersamp + k + 1;
            if k == 0; endi = 1; end
            set(0,'currentfigure',hlines)
            [starti,endi,yf,y] = xcorr_align(yf, Sf, S, endi, pntsperline);
            Yxcorr(line_num,1:min(pntsperline,length(y))) = y;
            if any(linexc == line_num)
              Y(line_num,1:min(pntsperline,length(y))) = y;
              align_ind(line_num,:) = align_ind(1,:);
            else
              set(0,'currentfigure',hlines)
              [dyflag,rise_temp,fall_temp] = find_rise(yf,line_num,top_thresh,bot_thresh);
              rtp = rise_temp;
              ftp = fall_temp;
              
              align_ind(line_num,:) = dyflag;
              if sum(dyflag) > 0
                if ~exist('align_base')
                  align_base = dyflag;
                end
                rise_diffv(line_num,:) = dyflag(1,:) - align_base;
                if abs(rise_diffv(line_num,1)) < min(pntsperline*0.05,25) % was using 10 in v5.2, changed to 20 in v5.21
                  starti = starti + rise_diffv(line_num,1);
                  endi = min(starti+pntsperline-1,length(S));
                  align_ind(line_num,:) = dyflag;
                else
                  align_ind(line_num,:) = align_ind(line_num-1,:);
                  %starti = starti + (align_ind(line_num,1) - align_ind(1,1));
                  %endi = min(starti+pntsperline-1,length(S));
                end
                
                if starti < 1
                  y = [S(1).*ones(1,-1*starti + 1),S(1:pntsperline + starti - 1)];
                  yf = [Sf(1).*ones(1,-1*starti + 1),Sf(1:pntsperline + starti - 1)];
                elseif rise_diffv(line_num,1) && k == 0
                  y = [y(starti:end),S(endi-(starti-1)+1:endi)];
                  yf = [yf(starti:end),Sf(endi-(starti-1)+1:endi)];
                else
                  y = S(starti:endi);
                  yf = Sf(starti:endi);
                end
                
                if length(y) < pntsperline
                  yf = [yf,yf(end).*ones(1,pntsperline - length(y))];
                  y = [y,y(end).*ones(1,pntsperline - length(y))];
                end
                ind(line_num,:) = [starti,endi];
                Y(line_num,1:min(pntsperline,length(y))) = y;
              else
                y = S(starti:endi);
                yf = Sf(starti:endi);
                starti = endi - round(line_buffer/2);
                endi = min(starti+pntsperline-1,length(S));
              end
            end
            
            if k == linespersamp-1
              endi = 1 + line_buffer;
            end
            if p
              set(0,'currentfigure',hlines)
              %plot(y,'r');
              if line_num > 1; plot(Y(line_num-1,:),'b'); end;
              hold on;
              plot(yf,'g'); hold off;
              title(sprintf('Sample %d, Line %d', l, line_num));
              legend(sprintf('Line %d',line_num-1),...%sprintf('Orignal %d',line_num),...
                sprintf('Aligned %d',line_num));
              set(0,'currentfigure',himg);
              imagesc(Y);
            end
          end
        end
      end
      for l = skipped_samples(:)'
        for k = 0:linespersamp-1
          line_num = (l-1)*linespersamp + k + 1;
          align_ind(line_num,:) = round(mean(align_ind(valid_lines(line_num:linespersamp:end),:),1));
          ind(line_num,:) = round(mean(ind(valid_lines(line_num:linespersamp:end),:),1));
        end
      end
      d = setfield(d,fldnm{f},Y);
    else
      d = setfield(d,fldnm{f},C);
    end
  end
  out_struct.Yxcorr = Yxcorr;
  out_struct.final_alind = align_ind;
  out_struct.ind = ind;
  out_struct.rise_diff = rise_diffv;
  out_struct.alskipped = alskipped;
  out_struct.thresh_perc = [top_perc,bottom_perc];
  varargout{1} = out_struct;
  
  for f = 1:length(fldnm)%
    C = getfield(din,fldnm{f});
    if isempty(strmatch(fldnm{f},{'line_datevec','dataset'}))
      C = C(:,1:pntspersamp);
      N = size(C,2);
      L = size(C,1);
      Y = zeros(L*linespersamp,pntsperline);
      for l = 1:L%eval_samp_num
        y = C(l,:);
        for k = 1:linespersamp
          i = ind((l-1)*linespersamp + k,1);
          yi = i:min(i+pntsperline-1,N);
          if i > 0
            if length(yi) < pntsperline
              yi = [yi,yi(end).*ones(1,pntsperline - length(yi))];
            end
            Y((l-1)*linespersamp + k,1:min(pntsperline,length(yi))) = y(yi);
          else
            Y((l-1)*linespersamp + k,1:min(pntsperline,length(yi))) = ...
              [y(1).*ones(1,-1*i + 1),y(1:pntsperline + i - 1)];
          end
        end
      end
      d = setfield(d,fldnm{f},Y);
    else
      d = setfield(d,fldnm{f},C);
    end
  end
  f = strmatch('line_datevec',fldnm);
  t = strmatch('Time',fldnm);
  C = getfield(din,fldnm{f});
  T = getfield(din,fldnm{t});
  dt = T(1,2)-T(1,1);
  if isnan(C)
    L = 1;
    N = 6;
  else
    L = size(C,1);
    N = size(C,2);
  end
  Y = zeros(L*linespersamp,N);
  for l = eval_samp_num
    %disp(l)
    if ~isempty(C) & ~isnan(C)
      y = C(l,:);
      yscaler = datenum(y);
    else
      dsetstr = getfield(d,'dataset');
      YYYY = str2num(dsetstr(1:4)); 
      MM = str2num(dsetstr(5:6)); 
      DD = str2num(dsetstr(7:8));
      yscaler = datenum([YYYY,MM,DD,00,00,00]);
    end
    tline = T(l,:);
    for k = 1:linespersamp
      i = ind((l-1)*linespersamp + k,1);
      %disp(k)
      if i > 0
        Y((l-1)*linespersamp + k,:) = datevec(yscaler + tline(i)/(3600*24));
      else
        Y((l-1)*linespersamp + k,:) = datevec(yscaler + (tline(1)+i*dt)/(3600*24));
      end
    end
  end
  d = setfield(d,fldnm{f},Y);
  ldv = d.line_datevec;
  d = rmfield(d,'dataset');
  if strmatch('a',keep_lines)
    keep_lines = 1:line_num;
  else
    d = structfun(@(x) (x(keep_lines,:)),d,'uniformoutput',false);
  end
  d = rmfield(d,'line_datevec');
  i1eq0 = find(ind(:,1) == 0);
  for q = 1:length(i1eq0)
    prob_ind = mod(i1eq0(q),10);
    if ~prob_ind; prob_ind = 10; end
    irelevant = prob_ind:10:size(ind,1);
    mem = ismember(irelevant,i1eq0(q));
    irelevant = irelevant(find(mem == 0));
    ind(i1eq0(q),1) = round(mean(ind(irelevant,1)));
  end
  iendeq0 = find(ind(:,2) == 0);
  for q = 1:length(iendeq0)
    prob_ind = mod(iendeq0(q),10);
    if ~prob_ind; prob_ind = 10; end
    irelevant = prob_ind:10:size(ind,1);
    mem = ismember(irelevant,iendeq0(q));
    irelevant = irelevant(find(mem == 0));
    ind(iendeq0(q),2) = round(mean(ind(irelevant,2)));
  end
  start_left = max(1,1 - min(ind(:,1)));
  l2lind_diff = ind(1:end-1,2)-ind(2:end,1);
  rm_right = max(abs(l2lind_diff(l2lind_diff(1:(keep_lines(end)-1)) < pntsperline*.25))); %max(abs(l2lind_diff(l2lind_diff(1:end-round(line_num*.05)) < pntsperline/2)));
  d = structfun(@(x) (x(:,start_left:pntsperline-rm_right)),d,'uniformoutput',false);
  d.dataset = din.dataset;
  d.line_datevec = ldv;
  
else
  disp('This dataset has already been reshaped');
  d = din;
end
close all
plot_heatmap(d,[],'scale',0,'thresh','auto','plot','lin'); %colormap(gray);

if svheat
  save_open_figures(d,[],[],'');
end

if svdat
  currentfldr = cd;
  Dpos = strfind(currentfldr,'Documents');

  home_dir = strcat(currentfldr(1:Dpos-1),'Documents\MADLab Data\Laser Ablation\',d.dataset,'\');
  allfiles = ls;
  files = cellstr(allfiles);
  
  for i = 1:length(fldnm)
    if isempty(strmatch(fldnm{i},{'dataset'}))
      csvname = strcat(home_dir,fldnm{i},'.csv');
      csvwrite(csvname,getfield(d,fldnm{i}));
    end
  end
end
%}

  function [starti,endi,yf, y] = xcorr_align(templine, S, So, endi, pntsperline);
    buffer = round(pntsperline/10);
    starti = endi-buffer;
    if starti <= 0; starti = 1; end
    endi = min(pntsperline+starti-1,length(S));
    matchline = S(starti:endi);
    unfilt_line = So(starti:endi);
    [R,lags] = xcorr(matchline(1:min(length(matchline),length(templine)))',...
      templine(1:min(length(matchline),length(templine)))');
    [m,i] = max(R);
    starti = starti + lags(i);
    endi = min(starti+pntsperline-1,length(S));
    if starti < 1
      yf = [matchline(1).*ones(1,-1*starti + 1),matchline(1:end + starti - 1)];
      y = [unfilt_line(1).*ones(1,-1*starti + 1),unfilt_line(1:end + starti - 1)];
      %starti = 1;
      %endi = min(starti+pntsperline-1,length(S));
    else
      yf = S(starti:endi);
      y = So(starti:endi);
    end
    if p
      plot(templine); hold on;
      plot(matchline,':m');
      plot(yf,':c');
    end
  end

end
