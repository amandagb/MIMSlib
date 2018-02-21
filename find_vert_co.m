function dout = find_vert_co(data,elem_ana_range,varargin);

%% Script: dout = find_vert_co(data,varargin);
% Description: Finds verticle cutoffs to separate background from tissue on
% left and right sides. Based on drift_correct_v3 function.
% Example:
% Required Functions: interp_data_elemrng, find_laoff_ind
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'Nma':  integer indicating  length of moving average filter. The longer
%   the filter, the greater the degree of smoothing on the mean sample
%   signal [DEFAULT = 30]
%   'lpfco': integer indicating LPF cutoff "frequency". [DEFUALT = 0.05]
%   'lpfo': FIR LPF order (number of coefficients [DEFUALT = 100]
%   'plot': binary indicating whether to generate plot (1) or not (0)
%   [DEFAULT = 1]
%   'lco': numeric indicating left cutoff index
%   'rco': numeric indicating right cutoff index
%   'laoff': matrix indicating the indices associated with laser off
% OUTPUTS ---------------------------------------------------------------
% dout:
%
%  Date           Author            E-mail                      Version
%  26 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if isstruct(data)
  fldnm = fieldnames(data);
  Ns = size(getfield(data,fldnm{1}),2);
else
  Ns = size(data,2);
end

if strmatch('laoff',PropertyNames)
  la0ind = PropertyVal{strmatch('laoff',PropertyNames)};
end

if strmatch('nma',PropertyNames)
  N_MA = PropertyVal{strmatch('nma',PropertyNames)};
else
  N_MA = round(Ns/10) - (mod(round(Ns/10),2)~=0);
  % Sets the moving average filter order to ~1/10 of the number of samples
  % and ensures that it is even since the mean signal is zero-padded on
  % either side by same number of zeros
end

if strmatch('lpfco',PropertyNames)
  lpf_co = PropertyVal{strmatch('lpfco',PropertyNames)};
else
  lpf_co = 0.05;
end

if strmatch('lpfo',PropertyNames)
  lpf_order = PropertyVal{strmatch('lpfo',PropertyNames)};
else
  lpf_order = round(Ns/5);
end

dout = [];
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
for f = elem_ana_range
  if isstruct(data)
    v = getfield(data,fldnm{f});
  else
    v = data;
  end
  
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    
    if ~isempty(find(isnan(v) == 1))
      nan_r = find(any(isnan(v),2)==1);
      v_non0 = v(find(any(isnan(v),2)==0),:);
      if isempty(v_non0)
        v_non0 = v(:,any(isnan(v),1));
      end
    else
      v_non0 = v;
    end
    
    % ---------- (1): Finds mean of sample along a line and applies a
    % moving average filter to the sample mean
    mu_S = mean(v_non0,1);
    h_MA = ones(1,N_MA)./N_MA;
    pad_mu_S = [ones(1,N_MA/2).*mu_S(1),mu_S,ones(1,N_MA/2-1).*mu_S(end)];
    mu_MA = conv(pad_mu_S,h_MA,'valid');
    
    % ---------- (2): Finds a threshold such that when applied to ?MA, the
    % resulting binary signal goes from 0 ? 1 ? 0
    [sort_mu_MA,sort_ind] = sort(mu_MA);
    
    i = 1;
    figure;
    while i
      thresh = sort_mu_MA(i);
      bv = zeros(1,length(mu_MA));
      bv(find(mu_MA > thresh)) = 1;
      d_bv = diff(bv);
      non0ind = find(d_bv ~= 0);
      plot([mu_S;mu_MA]','.');hold on; plot(bv.*max(mu_MA),'ro'); hold off;
      if length(non0ind) == 2 && non0ind(1)>Ns/10 && non0ind(2)<length(mu_MA)-Ns/10
        thresh = (max(mu_MA)-mu_MA(non0ind(1)))/3+mu_MA(non0ind(1));
        bv = zeros(1,length(mu_MA));
        bv(find(mu_MA > thresh)) = 1;
        d_bv = diff(bv);
        non0ind = find(d_bv ~= 0);
        i = 0;
      else
        i = i + 1;
      end
    end
    
    % ---------- (3): Calculate derivative of MAfiltered mean and then
    % smooth the result with a MA filter
    d_mu_MA = diff(mu_MA);
    MAd_mu_MA = conv([ones(1,N_MA/2).*d_mu_MA(1),d_mu_MA,ones(1,N_MA/2-1).*d_mu_MA(end)],h_MA,'valid');
    
    % ---------- (4): Calculate second derivative of MAfiltered mean then
    % low pass filter the result
    dd_mu_MA = diff(d_mu_MA);
    % LPFdd_mu_MA = conv([ones(1,N_MA/2).*dd_mu_MA(1),dd_mu_MA,ones(1,N_MA/2-1).*dd_mu_MA(end)],h_MA,'valid');
    dMAd_mu_MA = diff(MAd_mu_MA);
    LPF_coeff = fir1(lpf_order,lpf_co);
    LPFdMAd_mu_MA = filtfilt(LPF_coeff,1,dMAd_mu_MA);
    
    %-------------- Frequency Plots -------------------
    %     Fs = 1; %this is saying that time btwn samples is 1 second...this is not true but it really doesn't matter for now;
    %     L = length(d_mu_MA);
    %     NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    %     FTdmuMA = fft(d_mu_MA,NFFT)/L;
    %     FTMAdmuMA = fft(MAd_mu_MA,NFFT)/L;
    %     FTlpfMAdmuMA = fft(LPFdMAd_mu_MA,NFFT)/L;
    %     fvec = Fs/2*linspace(0,1,NFFT/2+1);
    %
    %     % Plot single-sided amplitude spectrum.
    %     fntsz = 14;
    %     figure('color','w');
    %     plot(fvec,2.*[abs(FTdmuMA(1:NFFT/2+1));abs(FTMAdmuMA(1:NFFT/2+1))],...
    %       'linewidth',2)
    %     %title('Fourier transform of line mean and median','fontsize',fntsz)
    %     xlabel('Frequency (Hz)','fontsize',fntsz);
    %     ylabel('Magnitude','fontsize',fntsz);
    %     %legend('FT\{mean\}','FT\{median\}','fontsize',fntsz);
    
    %-------------- Compile Output Variable -------------------
    mu_dout.mu_S = mu_S;
    mu_dout.d_mu_S = diff(mu_S);
    mu_dout.dd_mu_S = diff(diff(mu_S));
    mu_dout.mu_MA = mu_MA;
    mu_dout.d_mu_MA = d_mu_MA;
    mu_dout.dd_mu_MA = dd_mu_MA;
    mu_dout.MAd_mu_MA = MAd_mu_MA;
    mu_dout.dMAd_mu_MA = dMAd_mu_MA;
    mu_dout.LPFdMAd_mu_MA = LPFdMAd_mu_MA;
    mu_dout.params.LPF_cutoff = lpf_co;
    mu_dout.params.LPF_order = lpf_order;
    mu_dout.params.MA_order = N_MA;
    mu_dout.params.non0ind = non0ind;
    mu_dout.params.thresh = thresh;
    if ~exist('la0ind','var')
      la0ind = find_laoff_ind(v,fldnm{f});%mu_dout,fldnm{f});
    end
    
    if strmatch('plot',PropertyNames)
      gen_fig = PropertyVal{strmatch('plot',PropertyNames)};
    else
      gen_fig = 1;
    end
    
    % ---------- (5): Find the right and left cut off by finding the max
    % and min rate of change (max and min points of MAd_mu_MA) and then
    % search for max/min of second derivative within a given search range.
    % sample location of max/min will give l- and r-cutoffs
    srch_rng = 30;%non0ind(1)/2);%
    if la0ind(1) < length(mu_S)/2 %indicating LA off data comprises BEGINNING samples
      max_srch_start = max(la0ind(end)+1,round(non0ind(1) - srch_rng));
      max_srch_end = min(round(non0ind(1)+srch_rng),length(MAd_mu_MA));
      min_srch_start = max(1,round(non0ind(2)-srch_rng));
      min_srch_end = min(round(non0ind(2)+srch_rng),length(MAd_mu_MA));
      lco_srch_start = la0ind(end)+1;
      rco_srch_end = length(LPFdMAd_mu_MA);
    else
      max_srch_start = max(1,round(non0ind(1) - srch_rng));
      max_srch_end = min(round(non0ind(1)+srch_rng),length(MAd_mu_MA));
      min_srch_start = max(1,round(non0ind(2)-srch_rng));
      min_srch_end = min(round(non0ind(2)+srch_rng),la0ind(1)-1);
      lco_srch_start = 1;
      rco_srch_end = la0ind(1)-1;
    end
    [max_slope,imax_slope] = max(MAd_mu_MA(max_srch_start:max_srch_end));
    [min_slope,imin_slope] = min(MAd_mu_MA(min_srch_start:min_srch_end));
    imax_slope = imax_slope + max_srch_start;
    imin_slope = imin_slope + min_srch_start;
    dMAd_mu_MA = LPFdMAd_mu_MA;
    
    % Method for selecting right and left vertical cutoffs. SHOULD BE IMPROVED
    % IN SUBSEQUENT VERSIONS
    if strmatch('rco',PropertyNames)
      r_co = PropertyVal{strmatch('rco',PropertyNames)};
    else
      r_co = 0;
    end
    
    if strmatch('lco',PropertyNames)
      l_co = PropertyVal{strmatch('lco',PropertyNames)};
    else
      l_co = 0;
    end
    
    while ~(r_co*l_co)
      left_cutoff = find(dMAd_mu_MA == max(dMAd_mu_MA(lco_srch_start:imax_slope)));
      right_cutoff = find(dMAd_mu_MA == max(dMAd_mu_MA(imin_slope:rco_srch_end)));
      if mu_MA(left_cutoff) > thresh
        imax_slope = imax_slope - 1;
      else
        l_co = left_cutoff;
      end
      
      if mu_MA(right_cutoff) > thresh
        imin_slope = imin_slope + 1;
      else
        r_co = right_cutoff;
      end
    end
    
    mu_dout.la0ind = la0ind;
    mu_dout.l_co = l_co;
    mu_dout.r_co = r_co;
    
    %-------------- Plots -------------------
    if gen_fig
      fntsz = 14;
      figure('color','w','name',fldnm{f});
      [ax,h1,h2] = plotyy(1:length(mu_S),mu_S,1:length(mu_S)-1,MAd_mu_MA);
      set(h1,'linestyle','none','marker','.','displayname','Vertical Avg ($\mu_S%)');
      set(h2,'linewidth',2,'linestyle','-','color','c','marker','s',...
        'markersize',3,'markerfacecolor','c',...
        'displayname','Filtered $\frac{d\mu_{MA}}{dx}$ $\left($\frac{d\mu_{MA}}{dx}$\right)_{MA}$');
      hold(ax(1), 'on');
      set(ax(1),'fontsize',fntsz);
      ylabel(ax(1),'Intensity (cps)','fontsize',fntsz);
      plot(ax(1),1:length(mu_MA),mu_MA,'.g',...
        'displayname','Filtered $\mu_S$ ($\mu_{MA}$)');
      plot(ax(1),[0,length(mu_MA)],thresh.*[1,1],'--r','linewidth',2,...
        'displayname','Step Threshold');
      hold(ax(2), 'on');
      %plot(ax(2),1:length(mu_S)-1,d_mu_MA);
      plot(ax(2),[0,length(mu_MA)],[0,0],'k','linewidth',2,...
        'displayname','zero');
      scale_dMAd1 = min(max(MAd_mu_MA),abs(min(MAd_mu_MA)))/...
        max(max(dMAd_mu_MA),abs(min(dMAd_mu_MA)));
      scale_dMAd2 = 0;%abs(min(MAd_mu_MA)/min(dMAd_mu_MA));
      plot(ax(2),1:length(mu_S)-2,dMAd_mu_MA.*max(scale_dMAd1,scale_dMAd2)...
        ,'-m','linewidth',3,'marker','s','markersize',2,'markerfacecolor','m',...
        'displayname','LPF $\frac{d}{dx}\left(\left($\frac{d\mu_{MA}}{dx}$\right)_{MA}\right)$');
      set(ax(2),'xgrid','on','ygrid','on','fontsize',fntsz);
      xlabel(ax(2),'Sample','fontsize',fntsz);
      ylabel(ax(2),'Intensity (cps)','fontsize',fntsz);
      set(ax,'xlim',[0,length(mu_S)])
      yl1= get(ax(1),'ylim');
      yl2 = get(ax(2),'ylim');
      plot(ax(2),[l_co,l_co],yl2,'--r','linewidth',2);
      txpos = 50;
      text(l_co+2,yl1(2)-(yl1(2)-yl1(1))/txpos,num2str(l_co),'fontsize',fntsz);
      plot(ax(2),[r_co,r_co],yl2,'--r','linewidth',2);
      text(r_co+2,yl1(2)-(yl1(2)-yl1(1))/txpos,num2str(r_co),'fontsize',fntsz);
      plot(ax(2),[la0ind(1),la0ind(1)],yl2,'--r','linewidth',2);
      text(la0ind(1)+2,yl1(2)-(yl1(2)-yl1(1))/txpos,num2str(la0ind(1)),'fontsize',fntsz);
      %legend(ax(1),'interpreter','latex');
      title(sprintf('Background metrics for %s',fldnm{f}),'fontsize',fntsz);
      
      mu_dout.params.ax = ax;
      mu_dout.params.h1 = h1;
      mu_dout.params.h2 = h2;
    end
    
    dout = setfield(dout,fldnm{f},mu_dout);
  end
  dockf on all;
end
end