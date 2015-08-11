function [mask] = lxl_mask(data,elem_ana_range,varargin);

%% Script:  [mask] = lxl_mask(data,elem_ana_range,varargin);
% Description:  Defines line by line (lxl) mask. Uses simiar idea as that
% behind finding single and double pixel noise (as in fix_sngl_pxl).
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotope analyzed
% varargin - 'PropertyName','PropertyValue'
%   'winlen':  numeric specifying window length [DEFAULT = 5]
%   'p2psim':  numeric fraction indicating the degree of point to point
%   similarity requires, used to find spikes in data [DEFAULT = 0.7]
%   'lines': row matrix indicating specific lines user would like to
%   analyze [DEFAULT = 'all']
%   'plot': binary indicating whether to show plots (1) or not (0) [DEFAULT
%   = 0]
%   'npnts': number of consecutive high points to cut out [DEFAULT = 2]
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  4 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  5 Sept 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     2


mask = [];
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('winlen',PropertyNames)
  winlen = PropertyVal{strmatch('winlen',PropertyNames)};
else
  winlen = 5;
end

if strmatch('lin',PropertyNames) 
  lines = PropertyVal{strmatch('lin',PropertyNames)};
else
  lines = 1:size(getfield(data,fldnm{1}),1);
end

if strmatch('plot',PropertyNames) 
  p = PropertyVal{strmatch('plot',PropertyNames)};
  fntsz = 14;
else
  p = 0;
end

% bwareaopen: function to use to ensure contiguous points

for f = elem_ana_range
  if isstruct(data)
    d = getfield(data,fldnm{f});
  else
    d = data;
  end
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    row_col_ind = [];
    for l = lines
      lndat = d(l,:);
      lndiff = diff(lndat);
      pad_ln = [lndat(1).*ones(1,(winlen-1)/2),lndat,lndat(end).*ones(1,(winlen-1)/2)];
      pad_lndiff = [lndiff(1).*ones(1,(winlen+1)/2),lndiff,lndiff(end).*ones(1,(winlen-1)/2)];
      
      if p
        figure('color','w');
        plot([pad_ln;pad_lndiff]');
        hold on;
        legend('Line Intensity','Line Derivative');
        title(sprintf('%s Line %d',fldnm{f},l),'fontsize',fntsz)
      end
      
      ckind = [];
      fxind = [];
      for i = 1:length(pad_ln)-winlen
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
                if p; plot(g,pad_ln(g),'*r'); end
                pad_ln(g) = (pad_ln(g-1)+pad_ln(g+1))/2;
                if p; plot(g,pad_ln(g),'*g'); end
              elseif min(sum(der(1:npnts))*0.8,abs(der(end))) >= max(sum(der(1:npnts))*0.8,abs(der(end)))*p2psim ...
                  && npnts > 1
                fxind = [fxind,g-1,g];
                if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*r'); end
                y1 = pad_ln(g-npnts); y2 = pad_ln(g+1); x1 = g-npnts; x2 = g+1;
                m = (y2-y1)/(x2-x1);
                b = (x2*y1-x1*y2)/(x2-x1);
                pad_ln([g-1,g]) = m.*[g-1,g]+b;
                if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*g'); end
              end
            %% CASE 2 - checks if there's a single large spike and then a
            % large drop in intensity (single point). Searches for two
            % point sequence HIGH-LOW
            elseif der(end-1) > 0 && der(end) < 0
              if min(der(end-1),abs(der(end))) >= max(der(end-1),abs(der(end)))*p2psim
                fxind = [fxind,g];
                if p; plot(g,pad_ln(g),'*r'); end
                pad_ln(g) = (pad_ln(g-1)+pad_ln(g+1))/2;
                if p; plot(g,pad_ln(g),'*g'); end
              end
            %% CASE 3 - Checks if there's a large spike followed by two
            % neighboring sharp drops. Searching for three point sequence
            % HIGH-LOW-LOW
            elseif der(1) > 0 && sum(der(2:end) < 0) == length(der)-1
              if min(der(1),abs(sum(der(2:end))*0.8)) >= max(der(1),abs(sum(der(2:end))*0.8))*p2psim
                fxind = [fxind,g-1,g];
                if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*r'); end
                y1 = pad_ln(g-npnts); y2 = pad_ln(g+1); x1 = g-npnts; x2 = g+1;
                m = (y2-y1)/(x2-x1);
                b = (x2*y1-x1*y2)/(x2-x1);
                pad_ln([g-1,g]) = m.*[g-1,g]+b;
                if p; plot([g-1,g]',[pad_ln(g-1),pad_ln(g)]','*g'); end
              end
            end
          end
        end
      end
      d_corrected(l,:) = pad_ln(1+2:end-2);
      row_col_ind = [row_col_ind,[l.*(ones(1,length(fxind)));fxind]];
      if p; plot(pad_ln,'m'); end;
    end
    datafxed = setfield(datafxed,fldnm{f},d_corrected);
    fxedinfo = setfield(fxedinfo,fldnm{f},row_col_ind');
  elseif strmatch(fldnm{f},'Time')
    datafxed = setfield(datafxed,'Time',data.Time);
  elseif strmatch(fldnm{f},'line_datevec')
    datafxed = setfield(datafxed,'line_datevec',data.line_datevec);
  end
end
end