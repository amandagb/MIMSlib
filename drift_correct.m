function [dnorm,vco_out] = drift_correct(data,elem_ana_range,varargin);

%% Script: [dnorm,vco_out] = drift_correct(data,elem_ana_range,varargin);
% Description:
% Example: [dnorm,vco] =
% drift_correct(data,[],'plot',0,'stat',1,'freq',0,'heat',0,'norm','lpfmed'
% );
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'plot': binary indicating whether to generate plot (1) or not (0) in
%   find_vert_co function [DEFAULT = 1]
%   'heatmap': binary indicating whether to generate heatmap (1) or not (0)
%   [DEFUALT = 0]
%   'stats_plot': binary indicating whether to generate plot of image
%   statistics along vertical lines (1) or not (0) [DEFAULT = 1]
%   'freq_plot': binary indicating whether to generate plot of mean and
%   median fourier transform (1) or not (0) [DEFAULT = 0]
%   'vco': verticle cutoff structure generated by find_vert_co
%   'normalize_by': string indicating whether to subtract line mean, line
%   median, LPF line mean or LPF line median from each line. Inputs are
%   'mean','median','lpfmean','lpfmedian' [DEFUALT = 'mean']
% OUTPUTS ---------------------------------------------------------------
% dnorm: data structure which has same structure as that in data, but which
% is normalized by the
%
%  Date           Author            E-mail                      Version
%  24 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     Possible enchancements include (1) Give some automatic selection or
%     justification for values of Nma & fds, which will need to be adjusted
%     depending on the total number of samples, Ns; (2) Devise a better way
%     for selection r_co and l_co since as it stands, there are two
%     neighboring peaks in dd_mu_dsMA. Need to select the correct one.
%  25 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  26 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Added use of function find_laoff_ind to exclude LA indices
%  26 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Created find_vert_co function which returns verticle cutoff indices


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

dnorm = [];
stat = 1;
if strmatch('vco',PropertyNames)
  vco_out = PropertyVal{strmatch('vco',PropertyNames)};
else
  vco_out = [];
end

if ~isstruct(vco_out)
  vco_out = find_vert_co(data,elem_ana_range,varargin{:});
end

[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

for f = elem_ana_range
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    v = getfield(vco_out,fldnm{f});
    d = getfield(data,fldnm{f});
    bkgnd_ind = [1:v.l_co,v.r_co:size(d,2)];
    ex_ind = ismember(bkgnd_ind,v.la0ind);
    bkgnd_ind = bkgnd_ind(find(ex_ind == 0));
    
    if ~isempty(find(isnan(d) == 1))
      nan_r = find(any(isnan(d),2)==1);
      included_rows = find(any(isnan(d),2)==0);
      d_non0 = d(included_rows,:);
      if isempty(d_non0)
        d_non0 = d(:,any(isnan(d),1));
      end
    else
      d_non0 = d;
    end
    vbkgd = d_non0(:,bkgnd_ind);
    
    line_mean = mean(vbkgd,2);
    line_med = median(vbkgd,2);
    
    Fs = 1; %this is saying that time btwn samples is 1 second...this is not true but it really doesn't matter for now;
    L = length(line_mean);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    FTlnmean = fft(line_mean,NFFT)/L;
    FTlnmed = fft(line_med,NFFT)/L;
    fvec = Fs/2*linspace(0,1,NFFT/2+1);
    
    fplot = 0;
    if strmatch('freq',PropertyNames)
      fplot = PropertyVal{strmatch('freq',PropertyNames)};
    end
    
    if fplot
      % Plot single-sided amplitude spectrum.
      fntsz = 14;
      figure('color','w');
      plot(fvec,2.*[abs(FTlnmean(1:NFFT/2+1)),abs(FTlnmed(1:NFFT/2+1))])
      title('Fourier transform of line mean and median','fontsize',fntsz)
      xlabel('Frequency (Hz)','fontsize',fntsz);
      xlabel('Magnitude','fontsize',fntsz);
      legend('FT\{mean\}','FT\{median\}','fontsize',fntsz);
    end
    
    lpfco = 0.02; %LPF cutoff frequency in Hz
    lpfo = 10; % LPF filter order (suggested in MATLAB to be < length(signal)/3
    b = fir1(lpfo,2*pi*lpfco); %Gives FIR LPF coefficients
    lpfmean = filtfilt(b,1,line_mean);
    lpfmed = filtfilt(b,1,line_med);
    
    %  --------------------------------------------------------------------
    %   METHOD 3: LPF medians of each line -- developed 28 Aug 2011
    %  --------------------------------------------------------------------
    lpf_dat = [];
    if strmatch('norm',PropertyNames)
      norm_dat = PropertyVal{strmatch('norm',PropertyNames)};
      if strmatch('mean',norm_dat)
        norm_dat = line_mean;
        norm_str = 'Mean';
      elseif strmatch('med',norm_dat)
        norm_dat = line_med;
        norm_str = 'Median';
      elseif strmatch('lpf',norm_dat)
        if strmatch('lpfmean',norm_dat)
          norm_dat = lpfmean;
          norm_str = 'LPF Mean';
        elseif strmatch('lpfmed',norm_dat)
          norm_dat = lpfmed;
          norm_str = 'LPF Median';
        end
      end
    else
      norm_dat = lpfmed;
      norm_str = 'LPF Median';
    end
    
    %dnon0_norm = d_non0./repmat(LPFlndat,1,size(d_non0,2));  <== commented
    %out 29 Aug because of testing on 120710Brain 'Mo95' data. Multiple
    %datum had intensity values of 0 resulting in huge numbers in the
    %division. Method instead subtracts mean from all numbers in a given
    %line, then ensures there are no 0 values.
    %
    % QUESTION REMAINS AS TO WHETHER line_mean, line_median, LPFlnmed, or
    % LPFlnmean should be subtracted.
    dnon0_norm = d_non0 - repmat(norm_dat,1,size(d_non0,2));
    if min(dnon0_norm(:)) < 0
      dnon0_norm = dnon0_norm + abs(min(dnon0_norm(:)));
    end
    d_norm = nan(size(d));
    d_norm(included_rows,:) = dnon0_norm;
    if strmatch('heat',PropertyNames)
      heatmap = PropertyVal{strmatch('heat',PropertyNames)};
      if heatmap
        plot_heatmap(d_norm,sprintf('%s Normalized by %s',fldnm{f},norm_str),...
          'thresh','auto','plot','lin');
      end
    end
    dnorm = setfield(dnorm,fldnm{f},d_norm);
    
    if strmatch('stat',PropertyNames)
      stat = PropertyVal{strmatch('stat',PropertyNames)};
    end
    
    if stat
      fntsz = 14;
      figure('color','w'); plot(included_rows,[line_mean,line_med,lpfmean,lpfmed,...
        mean(d_norm(included_rows,bkgnd_ind),2),...
        median(d_norm(included_rows,bkgnd_ind),2)]);%,'-*');
      set(gca,'ygrid','on','fontsize',fntsz);
      xlabel('Line Number','fontsize',fntsz);
      ylabel('Intensity (cps)','fontsize',fntsz);
      legend('mean','median','LPFmedian','LPF median','d_{norm} mean','d_{norm} median',...
        'location','southeast');
      title(sprintf('Background metrics for %s -- Normalized by %s',fldnm{f},norm_str),...
        'fontsize',fntsz);
    end
    
  elseif strmatch(fldnm{f},'Time')
    dnorm = setfield(dnorm,'Time',data.Time);
  elseif strmatch(fldnm{f},'line_datevec')
    dnorm = setfield(dnorm,'line_datevec',data.line_datevec);
  end
  
  %  -----------------------------------------------------------------------
  %   METHOD 2: attempt to threshold very large values -- developed 28 Aug
  %   2011
  %  -----------------------------------------------------------------------
  %   [t,id] = auto_thresh(vbkgd,'auto',[]);
  %   [rsort,rid] = sort(id{2}(:,1));
  %   csort = id{2}(rid,2);
  %   unir = unique(rsort);
  %   for r = unir'
  %     disp(r);
  %     c = csort(rsort == r);
  %     lc = length(c);
  %     d_c = diff(c);
  %     non_conseq = [0,find(d_c ~= 1)',lc];
  %     %figure; plot(vbkgd(r,:));
  %     for n = 1:length(non_conseq)-1
  %       if c(non_conseq(n)+1) == 1
  %         c1 = 1;
  %         i1 = t(2);
  %       else
  %         c1 = c(non_conseq(n)+1) - 1;
  %         i1 = vbkgd(r,c1);
  %       end
  %
  %       if c(non_conseq(n+1))+1 > size(vbkgd,2)
  %         c2 = size(vbkgd,2);
  %         i2 = t(2);
  %       else
  %         c2 = c(non_conseq(n+1)) + 1;
  %         i2 = vbkgd(r,c2);
  %       end
  %
  %       m = (i1-i2)/(c1-c2);
  %       b = (c1*i2-c2*i1)/(c1-c2);
  %       vbkgd(r,c(non_conseq(n)+1:non_conseq(n+1))) = m.*c(non_conseq(n)+1:non_conseq(n+1)) + b;
  %     end
  %     %hold on; plot(vbkgd(r,:),'r');
  %   end
  %
  %   line_meanlinfit = mean(vbkgd,2);
  %   line_medlinfit = median(vbkgd,2);
  %   figure; plot(included_rows,[line_mean,line_meanlinfit,line_med,line_medlinfit]);
  %   legend('mean','Linear fit mean','median','Linear fit median','location','northwest');
  
  %  -----------------------------------------------------------------------
  %   METHOD 1: attempt to LPF values -- developed 27 Aug 2011
  %  -----------------------------------------------------------------------
  %   ln = 143-length(nan_r);
  %   lndat = vbkgd(ln,:);
  %   Fs = 1; %this is saying that time btwn samples is 1 second...this is not true but it really doesn't matter for now;
  %   L = length(bkgnd_ind);
  %   NFFT = 2^nextpow2(L); % Next power of 2 from length of y
  %   FTlndat = fft(lndat,NFFT)/L;
  %   f = Fs/2*linspace(0,1,NFFT/2+1);
  %
  %   % Plot single-sided amplitude spectrum.
  %   figure;
  %   plot(f,2*abs(FTlndat(1:NFFT/2+1)))
  %
  %   lpfco = 0.001; %in Hz
  %   lpfo = 10;
  %   b = fir1(lpfo,2*pi*lpfco);
  %   LPFlndat = filtfilt(b,1,lndat);
  %   LPFd_non0 = filtfilt(b,1,d_non0(:,[bkgnd_ind])');
  %   LPFd_non0 = LPFd_non0';
  %
  %   figure; plot([lndat;LPFd_non0(ln,:)]'); dockf on all
  %   LPFlnavg = mean(LPFd_non0,2);
  %   LPFlnmed = median(LPFd_non0,2);
  %   figure; plot(find(any(isnan(d),2)==0),[line_mean,LPFlnavg,line_med,LPFlnmed]);
  %   legend('mean','LPF mean','median','LPF median','location','northwest');
  
end
dockf on all
end