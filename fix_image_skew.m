function [d,varargout] = fix_image_skew(din,varargin)
%% Script: d = fix_image_skew(dstand,varargin)
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
%   'lines': integer indicating number of lines per sample [DEFAULT =
%   find_linespersamp function]
%   'temp_el': string indicate the template element to be used. The string
%   must uniquely identify the element and must have the first letter
%   capitalized. [DEFAULT = 'C_']
%   'plot': binary indicating whether to generate debugging plots (1) or
%   not (0) [DEFAULT = 0]
%   'save_heat': binary indicating whether to save heatmaps (1) or not (0)
%   [DEFAULT = 0]
%   'save_data': binary indicating whether to save restructured data (1) or
%   not (0) [DEFAULT = 0]
%   'eval_samp_num': evaluates specific sample numbers [DEFAULT = evaluate
%   all sample numbers]
%   'temp_lines': determine which lines to use for the template
%   'keep_lines': vector containing lines that should be kept of final
%   image [DEFAULT = 'all']
%   'nalign': indicates the number of points that are allowed to vary
%   between the rises of subsequent lines [DEFAULT = 20]
%   'top_perc': integer indicating the percentage of data that should be
%   above the upper threshold. [DEFAULT = automated selection*--still in
%   debugging stages as of v5.4, 23 May 2012--use 20120320_OESCugrid40umps
%   for validation. The automated selection is done in lines 226-248]
%   'bottom_perc': integer indicating the percentage of data that should be
%   below the lower threshold. [DEFAULT = automated selection*--still in
%   debugging stages as of v5.4, 23 May 2012--use 20120320_OESCugrid40umps
%   for validation. The automated selection is done in lines 202-224]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%   9 Nov  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%       CODE BASED ON reshape_OEScorr v5.6 but it modified to fix data that
%       has 1 "sample per line" and is simply wrapped data. The data should
%       be surrounded on either side by vallies with plateu in between.

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
himg = figure('position',[53   211   560   420]);%('position',[ 25   475   560   420]);%390   573   560   420]);%[403   246   560   420]);%
hlines = figure('position',[ 635   211   560   420]);%('position',[ 605   475   560   420]);%390    60   560   420]);%[403  -268   560   420]);%

if din.Time(1,1) == din.Time(2,1) % This condition indicates that the dataset has not previously been reshaped
  %% Generate Pattern using template element
  for f = interp_data_elemrng(din,temp_el)
    C = getfield(din,fldnm{f});
    if isempty(strmatch(fldnm{f},{'Time','line_datevec','dataset'}))
      C = C(:,1:pntspersamp);
      N = size(C,2); % number of samples
      L = size(C,1);  % number of lines
      Cf = medfilt1(C,5,L,2); % median filtering for large spike elimination by taking the median of 5 points
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
      line_buffer = 0;%round(pntsperline/10); % used as a line buffer
      starti = 1;
      endi = pntspersamp;
      %pntsperline = endi;
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
            %Yxcorr(line_num,1:min(pntsperline,length(y))) = y;
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
                starti = starti + rise_diffv(line_num,1);%1;%
                endi = min(starti+pntsperline-1,length(S));%pntspersamp;%
                align_ind(line_num,:) = dyflag;
                
                if starti < 1
                  y = [S(1).*ones(1,-1*starti + 1),S(1:pntsperline + starti - 1)];
                  yf = [Sf(1).*ones(1,-1*starti + 1),Sf(1:pntsperline + starti - 1)];
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
          align_ind(line_num,:) = round(mean(align_ind(valid_lines,:),1));
          ind(line_num,:) = round(mean(ind(valid_lines,:),1));
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
          %disp(sprintf('Line: %d; Field %s',l,fldnm{f}));
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
  L = size(C,1);
  N = size(C,2);
  Y = zeros(L*linespersamp,N);
  for l = eval_samp_num
    y = C(l,:);
    tline = T(l,:);
    yscaler = datenum(y);
    for k = 1:linespersamp
      i = ind((l-1)*linespersamp + k,1);
      if i > 0
        Y((l-1)*linespersamp + k,:) = datevec(yscaler + tline(i)/(3600*24));
      else
        Y((l-1)*linespersamp + k,:) = datevec(yscaler + (tline(1)+i*dt)/(3600*24));
      end
    end
  end
  d = setfield(d,fldnm{f},Y);
  
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
  home_dir = strcat('D:\My Documents\MADLab Data\Laser Ablation\',d.dataset,'\');
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

  function [starti,endi,corrlinef,corrline] = xcorr_align(templine, S, So, endi, pntsperline);
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
      corrlinef = [mean(matchline(1:5)).*ones(1,-1*starti + 1),matchline(1:end)];
      corrline = [mean(unfilt_line(1:5)).*ones(1,-1*starti + 1),unfilt_line(1:end)];
      %starti = 1;
      %endi = min(starti+pntsperline-1,length(S));
    else
      corrlinef = S(starti:end);
      corrline = So(starti:end);
    end
    if p
      plot(templine); hold on;
      plot(matchline,':m');
      plot(corrlinef,':c');
    end
  end

end
