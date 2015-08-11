function d = reshape_OESstructure(dstand,varargin)
%% Script: d = reshape_OESstructure(dstand,varargin)
% Description: Takes the OES standard text file with structure
% LABELline#,MM/DD/YYYY,
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'linespersample':
%   'temp_el': string indicate the template element to be used. The string
%   must uniquely identify the element and must have the first letter
%   capitalized. [DEFAULT = 'C_']
%   'rise_template_indices': matrix containing the index range containing
%   the rise template data [DEFAULT = 20:45]
%   'fall_template_indices': matrix containing the index range containing
%   the fall template data [DEFAULT = 220:250]
%   'plot': binary indicating whether to generate debugging plots (1) or
%   not (0) [DEFAULT = 0]
%   'save_heat': binary indicating whether to save heatmaps (1) or not (0)
%   [DEFAULT = 0]
%   'save_data': binary indicating whether to save restructured data (1) or
%   not (0) [DEFAULT = 0]
%   'eval_samp_num': evaluates specific sample numbers [DEFAULT = evaluate
%   all sample numbers]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  27 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  29 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Added ability to restructure the line_datevec and Time data
%  01 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Attempting to remove convolution from best template matching
%     computation
%  06 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Due to inability to analyze 20120228_OESpancreas dataset, a new
%     technique for alignment was developed using line gradient alignment
%  11 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     5

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

fldnm = fieldnames(dstand);

if strmatch('temp_el',PropertyNames)
  temp_el = PropertyVal{strmatch('temp_el',PropertyNames)};
else
  temp_el = 'C_';
end

if strmatch('lines',PropertyNames)
  linespersamp = PropertyVal{strmatch('lines',PropertyNames)};
else
  lindat = getfield(dstand,fldnm{interp_data_elemrng(dstand,temp_el)});
  lindat = lindat(1,:);
  T = mean(lindat)*0.9;
  indgtT = find(lindat > T);
  compind = diff(indgtT);
  linespersamp = length(find(compind > 30)) + 1;
end

if strmatch('rise',PropertyNames)
  rise_tmp_rng = PropertyVal{strmatch('rise',PropertyNames)};
else
  rise_tmp_rng = 20:45;
end

if strmatch('fall',PropertyNames)
  fall_tmp_rng = PropertyVal{strmatch('fall',PropertyNames)};
else
  fall_tmp_rng = 220:250;
end

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end

if strmatch('save_heat',PropertyNames)
  s = PropertyVal{strmatch('save_heat',PropertyNames)};
else
  s = 0;
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 0;
end

if strmatch('eval_sam',PropertyNames)
  eval_samp_num = PropertyVal{strmatch('eval_sam',PropertyNames)};
else
  eval_samp_num = 1:size(dstand.Time,1);
end

if strmatch('temp_lines',PropertyNames)
  temp_lines = PropertyVal{strmatch('temp_lines',PropertyNames)};
else
  temp_lines = 1:size(dstand.Time,1);
end

pntspersamp = size(dstand.Time,2) - 1;
if strmatch('pnts',PropertyNames)
  pntsperline = PropertyVal{strmatch('pnts',PropertyNames)};
else
  pntsperline = round(pntspersamp/linespersamp);
end
d = [];
flag_lines = [];

if dstand.Time(1,1) == dstand.Time(2,1)
  %% Generate Pattern
  for f = interp_data_elemrng(dstand,temp_el)
    C = getfield(dstand,fldnm{f});
    if isempty(strmatch(fldnm{f},{'Time','line_datevec','dataset'}))
      C = C(:,1:pntspersamp);
      N = size(C,2);
      L = size(C,1);
      Y = zeros(L*linespersamp,pntsperline);
      ind = zeros(L*linespersamp,2);
      rise_temp = mean(C(temp_lines,rise_tmp_rng));
      fall_temp = mean(C(temp_lines,fall_tmp_rng));
      for l = eval_samp_num
        y = [C(l,:),C(l,end).*ones(1,length(rise_temp)-1)];
        rise_temp_diff = zeros(N,1);
        for k = 1:N
          rise_temp_diff(k) = sum(rise_temp - y(k:k+length(rise_temp)-1));
        end
        rise_sign = diff(sign(rise_temp_diff));
        rise_starts = sort([find(rise_temp_diff == 0);...
          find(rise_sign == -2)]);
        y = y(1:N);
        
        if length(rise_starts) > linespersamp
          rise_search_start = [];
          risemaxthresh = 0.8;
          while length(rise_search_start) < linespersamp
            risemaxthresh = risemaxthresh - 0.05;
            rise_pos = find(rise_temp_diff >= max(rise_temp_diff)*risemaxthresh);
            rise_diffgtT = diff(rise_pos);
            rise_search_start = rise_pos(find(rise_diffgtT > 1));
          end
          rise_search_start_mat = repmat(rise_search_start',length(rise_starts),1);
          rise_starts_mat = repmat(rise_starts,1,length(rise_search_start));
          rise_row_diff = rise_starts_mat-rise_search_start_mat;
          rise_reprows = sum(diff(sign(rise_row_diff),1),2);
          r = 1;
          rs = [];
          for i = [find(rise_reprows)',length(rise_starts)]
            c = 1 + (sum( rise_reprows(1:min(i,length(rise_starts)-1) ) )...
              - rise_reprows(min(i,length(rise_starts)-1)))/2 ...
              + 0^(i < length(rise_starts))...
              - 0^(rise_reprows(end))*(i == length(rise_starts));
            disp(c)
            [v,k] = min(abs(rise_row_diff(1:i,c)));
            rs(r) = rise_starts(k);
            r = r + 1;
          end
          rise_starts = rs';
        end
        %rise_starts = rise_starts -length(rise_temp)+1;
        rise_starts(rise_starts < 0) = rise_starts(rise_starts < 0) + length(rise_temp) - 1;
        
        y = [C(l,:),C(l,end).*ones(1,length(fall_temp)-1)];
        fall_temp_diff = zeros(N,1);
        for k = 1:N
          fall_temp_diff(k) = sum(fall_temp - y(k:k+length(fall_temp)-1));
        end
        fall_sign = diff(sign(fall_temp_diff));
        fall_starts = sort([find(fall_temp_diff == 0);...
          find(fall_sign == 2)]);
        y = y(1:N);
        if length(fall_starts) > linespersamp
          fallmaxthresh = 0.7;
          fall_pos = find(fall_temp_diff >= max(fall_temp_diff)*fallmaxthresh);% 0.8); -- this number was changed on 2/28 TO 0.6 FROM 0.8
          fall_diffgtT = diff(fall_pos);
          fall_search_start = [fall_pos(find(fall_diffgtT > 1));fall_pos(end)];
          
          fall_search_start_mat = repmat(fall_search_start',length(fall_starts),1);
          fall_starts_mat = repmat(fall_starts,1,length(fall_search_start));
          fall_row_diff = fall_search_start_mat-fall_starts_mat;
          fall_reprows = sum(diff(sign(fall_row_diff),1),2);
          r = 1;
          fs = [];
          for i = [find(fall_reprows)',length(fall_starts)]
            c = 1 + (-sum( fall_reprows(1:min(i,length(fall_starts)-1) ) )...
              + fall_reprows(min(i,length(fall_starts)-1)))/2 ...
              + 0^(i < length(fall_starts))...
              - 0^(-fall_reprows(end))*(i == length(fall_starts)) + 1;
            [v,k] = min(abs(fall_row_diff(1:i,c)));
            fs(r) = fall_starts(k);
            r = r + 1;
          end
          fall_starts = fs';
        end
        %fall_starts = fall_starts -length(fall_temp)+1;
        
        %% Plotting routines
        if p
          Ivec = zeros(1,length(y));
          Ivec(rise_starts) = 1;
          rise_rep_temp = conv(Ivec,rise_temp);
          
          Ivec = zeros(1,length(y));
          Ivec(fall_starts) = 1;
          fall_rep_temp = conv(Ivec,fall_temp);
          %
          %           figure; plot(y./max(y),'r'); hold on;
          %           plot(rise_temp_diff./max(rise_temp_diff),'b')
          %           plot(rise_sign,'c')
          %           axis tight;
          %           figure; plot(y./max(y),'r'); hold on;
          %           plot(fall_temp_diff./max(fall_temp_diff),'b')
          %           plot(fall_sign,'c')
          %           axis tight; fallx = xlim;
          %           plot(fallx,[fallmaxthresh,fallmaxthresh],':k');
          figure; plot(y); hold on;
          plot(rise_rep_temp,'.g'); axis tight
          plot(rise_starts,rise_temp(1).*ones(1,length(rise_starts)),'oc');
          title(sprintf('Sample %d',l));
          plot(fall_rep_temp,'.r'); axis tight
          plot(fall_starts,fall_temp(1).*ones(1,length(fall_starts)),'om');
          title(sprintf('Sample %d',l));
        end
        %=================================================================
        buffer_pnts = 10;
        line_ends = fall_starts + length(fall_temp) + 10;
        disp(l)
        if length(rise_starts) == length(fall_starts) &&...
            length(rise_starts) == linespersamp
          disp([rise_starts,fall_starts])
          for k = 1:length(rise_starts)
            rng = rise_starts(k) - buffer_pnts : line_ends(k);
            rng = rng(rng > 0);
            yi = rng(1):min(rng(1)+pntsperline-1,N);
            ind((l-1)*linespersamp + k,:) = [yi(1),yi(end)];
            Y((l-1)*linespersamp + k,1:min(pntsperline,length(yi))) = y(yi);
          end
          if length(rise_starts) < linespersamp
            flag_lines = [flag_lines,l];
          end
        else
          flag_lines = [flag_lines,l];
        end
        dockf on all
      end
      d = setfield(d,fldnm{f},Y);
    else
      d = setfield(d,fldnm{f},C);
    end
  end
  
  if ~isempty(flag_lines)
    avg_lines_ind = eval_samp_num(find(ismember(eval_samp_num,flag_lines) == 0));
    for f = 1:length(flag_lines)
      for k = 1:linespersamp
        valid_rise_ind = (avg_lines_ind-1).*linespersamp + k;
        avg_rise_ind = round(mean(ind(valid_rise_ind,1)));
        valid_fall_ind = (avg_lines_ind-1).*linespersamp + k;
        avg_fall_ind = round(mean(ind(valid_fall_ind,2)));
        ind((flag_lines(f)-1).*linespersamp + k,:) = [avg_rise_ind, avg_fall_ind];
      end
    end
  end
  
  if strmatch(dstand.dataset,'20120228_OESpancreasPOS')
    ind([84,85],:) = ind([84,85],:) - 16;
  end
  
  for f = 1:length(fldnm)%
    C = getfield(dstand,fldnm{f});
    if isempty(strmatch(fldnm{f},{'line_datevec','dataset'}))
      C = C(:,1:pntspersamp);
      N = size(C,2);
      L = size(C,1);
      Y = zeros(L*linespersamp,pntsperline);
      for l = eval_samp_num
        y = C(l,:);
        for k = 1:linespersamp
          i = ind((l-1)*linespersamp + k,1);
          yi = i:min(i+pntsperline-1,N);
          if i
            Y((l-1)*linespersamp + k,1:min(pntsperline,length(yi))) = y(yi);
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
  C = getfield(dstand,fldnm{f});
  T = getfield(dstand,fldnm{t});
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
      if i
        Y((l-1)*linespersamp + k,:) = datevec(yscaler + tline(i)/(3600*24));
      end
    end
  end
  d = setfield(d,fldnm{f},Y);
else
  disp('This dataset has already been reshaped');
  d = dstand;
end

close all
plot_heatmap(d,[],'scale',0,'thresh','auto','plot','lin'); %colormap(gray);

if s
  h = get(0,'children');
  set(h,'PaperPositionMode','auto');
  fign = get(h,'Name');
  for i=h'
    %set(h,'colormap',gray);
    imgnm = sprintf('D:\\My Documents\\MADLab Data\\Laser Ablation\\%s\\Images\\%s %s',...
      d.dataset, d.dataset, fign{i}(1:6));
    saveas(h(i), imgnm, 'jpg');
  end
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
end
