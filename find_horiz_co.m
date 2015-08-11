function dout = find_horiz_co(data,elem_ana_range,varargin);

%% Script: dout = find_horiz_co(data,elem_ana_range,varargin);
% Description: Finds the upper and lower horizontal cutoff lines which
% separate background from tissue lines (horizontal lines). Tested on
% 120710Brain data with Nma = 6, lpf_co = 0.2, lpf_order = 6, srch_rng =
% 10. Works well for all except Ni & some problems with Cr.
% Example: hco = find_horiz_co(data,[],'dnorm',dnorm,'vco',vco);
% Required Functions: interp_data_elemrng, drift_correct
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
% OUTPUTS ---------------------------------------------------------------
% dout: data structure which has same first level structure as that in
% data, but with additional fields. Most important/useful fields are
% dout.ELEMENT{.up_co} and {.low_co} which tell the upper and lower
% horizontal cutoffs of the data.
%
%  Date           Author            E-mail                      Version
%  29 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('nma',PropertyNames)
  N_MA = PropertyVal{strmatch('nma',PropertyNames)};
else
  N_MA = 6;
end

if strmatch('lpfco',PropertyNames)
  lpf_co = PropertyVal{strmatch('lpfco',PropertyNames)};
else
  lpf_co = 0.2;
end

if strmatch('lpfo',PropertyNames)
  lpf_order = PropertyVal{strmatch('lpfo',PropertyNames)};
else
  lpf_order = 6;
end

if strmatch('vco',PropertyNames)
  vco_out = PropertyVal{strmatch('vco',PropertyNames)};
else
  vco_out = [];
end

if strmatch('dnorm',PropertyNames)
  dnorm = PropertyVal{strmatch('dnorm',PropertyNames)};
else
  dnorm = [];
end

if isempty(dnorm) | isempty(vco_out)
  [dnorm,vco_out] = drift_correct(data,elem_ana_range,varargin{:});
end

[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
dout = [];
for f = elem_ana_range
  v = getfield(dnorm,fldnm{f});
  
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
    
    mu_S = mean(v_non0,2)';
    h_MA = ones(1,N_MA)./N_MA;
    pad_mu_S = [ones(1,N_MA/2).*mu_S(1),mu_S,ones(1,N_MA/2-1).*mu_S(end)];
    mu_MA = conv(pad_mu_S,h_MA,'valid');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     lndat = mu_S;
    %     Fs = 1; %this is saying that time btwn samples is 1 second...this is not true but it really doesn't matter for now;
    %     L = length(mu_S);
    %     NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    %     FTlndat = fft(lndat,NFFT)/L;
    %     fvec = Fs/2*linspace(0,1,NFFT/2+1);
    %
    %     % Plot single-sided amplitude spectrum.
    %     figure;
    %     plot(fvec,2*abs(FTlndat(1:NFFT/2+1)))
    %
    %     lpfco = 0.0005; %in Hz
    %     lpfo = 10;
    %     b = fir1(lpfo,2*pi*lpfco);
    %     LPFlndat = filtfilt(b,1,lndat);
    %
    %     figure; plot([lndat;LPFlndat]');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [sort_mu_MA,sort_ind] = sort(mu_MA);
    
    i = 1;
    while i
      thresh = sort_mu_MA(i);
      bv = zeros(1,length(mu_MA));
      bv(find(mu_MA > thresh)) = 1;
      d_bv = diff(bv);
      non0ind = find(d_bv ~= 0);
      if length(non0ind) == 2 && non0ind(1)>10 && ...
          non0ind(1) < length(mu_S)/2 && non0ind(2)<length(mu_MA)-10
        thresh = thresh;%(max(mu_MA)-mu_MA(non0ind(1)))/3+mu_MA(non0ind(1));
        bv = zeros(1,length(mu_MA));
        bv(find(mu_MA > thresh)) = 1;
        d_bv = diff(bv);
        non0ind = find(d_bv ~= 0);
        i = 0;
      elseif i == length(sort_mu_MA)
        i = [];
        thresh = (max(mu_MA)-min(mu_MA))/3+min(mu_MA);
        non0ind = round([1+length(mu_MA)/4,length(mu_MA)-length(mu_MA)/4]);
        bv(:) = 1;
        bv([1:non0ind(1),non0ind(2):length(bv)])=0;
      else
        i = i + 1;
      end
    end
    
    [t,id] = auto_thresh(mu_S,'auto',10);
    if ~isempty(find(bv == 0 & mu_S > thresh))
      c = find(bv == 0 & mu_S > thresh);
      d_c = diff(c);
      lc = length(c);
      non_conseq = [0,find(d_c ~= 1),lc];
      %figure; plot(vbkgd(r,:));
      for n = 1:length(non_conseq)-1
        if c(non_conseq(n)+1) == 1
          c1 = c(non_conseq(n+1)) + 1;
          i1 = mu_S(c1);
        else
          c1 = c(non_conseq(n)+1) - 1;
          i1 = mu_S(c1);
        end
        
        if c(non_conseq(n+1))+1 > length(mu_S)
          c2 = c1;
          i2 = mu_S(c2);
        else
          c2 = c(non_conseq(n+1)) - 1;
          i2 = mu_S(c2);
        end
        
        m = (i1-i2)/(c1-c2);
        b = (c1*i2-c2*i1)/(c1-c2);
        
        if isnan(m)
          mu_S(c(non_conseq(n)+1):c(non_conseq(n+1)))= i2;
        else
          mu_S(c(non_conseq(n)+1):c(non_conseq(n+1)))= m.*[c(non_conseq(n)+1):c(non_conseq(n+1))] + b;
        end
      end
        h_MA = ones(1,N_MA)./N_MA;
        pad_mu_S = [ones(1,N_MA/2).*mu_S(1),mu_S,ones(1,N_MA/2-1).*mu_S(end)];
        mu_MA = conv(pad_mu_S,h_MA,'valid');
    end
    
    d_mu_MA = diff(mu_MA);
    MAd_mu_MA = conv([ones(1,N_MA/2).*d_mu_MA(1),d_mu_MA,ones(1,N_MA/2-1).*d_mu_MA(end)],h_MA,'valid');
    
    dd_mu_MA = diff(d_mu_MA);
    % LPFdd_mu_MA = conv([ones(1,N_MA/2).*dd_mu_MA(1),dd_mu_MA,ones(1,N_MA/2-1).*dd_mu_MA(end)],h_MA,'valid');
    dMAd_mu_MA = diff(MAd_mu_MA);
    LPF_coeff = fir1(lpf_order,lpf_co);
    LPFdMAd_mu_MA = filtfilt(LPF_coeff,1,dMAd_mu_MA);
    
    mu_dout.mu_L = mu_S;
    mu_dout.d_mu_L = diff(mu_S);
    mu_dout.dd_mu_L = diff(diff(mu_S));
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
    
    if strmatch('plot',PropertyNames)
      gen_fig = PropertyVal{strmatch('plot',PropertyNames)};
    else
      gen_fig = 1;
    end
    
    srch_rng = 10;%non0ind(1)/2);%
    max_srch_start = max(1,round(non0ind(1) - srch_rng));
    max_srch_end = min(round(non0ind(1)+srch_rng),length(MAd_mu_MA));
    min_srch_start = max(1,round(non0ind(2)-srch_rng));
    min_srch_end = min(round(non0ind(2)+srch_rng),length(MAd_mu_MA));
    upco_srch_start = 1;
    lowco_srch_end = length(LPFdMAd_mu_MA);
    
    [max_slope,imax_slope] = max(MAd_mu_MA(max_srch_start:max_srch_end));
    [min_slope,imin_slope] = min(MAd_mu_MA(min_srch_start:min_srch_end));
    imax_slope = imax_slope + max_srch_start-1;
    imin_slope = imin_slope + min_srch_start-1;
    dMAd_mu_MA = LPFdMAd_mu_MA;
    
    % Method for selecting lower and upper vertical cutoffs. SHOULD BE IMPROVED
    % IN SUBSEQUENT VERSIONS
    up_co = 0; low_co = 0;
    while ~(low_co*up_co)
      upper_cutoff = find(dMAd_mu_MA == max(dMAd_mu_MA(upco_srch_start:imax_slope)));
      lower_cutoff = find(dMAd_mu_MA == max(dMAd_mu_MA(imin_slope:lowco_srch_end)));
      if mu_MA(upper_cutoff) > thresh
        imax_slope = imax_slope - 1;
      else
        up_co = upper_cutoff;
      end
      
      if mu_MA(lower_cutoff) > thresh
        imin_slope = imin_slope + 1;
      else
        low_co = lower_cutoff;
      end
    end

    mu_dout.up_co = up_co;
    mu_dout.low_co = low_co;
    
    if gen_fig
      fntsz = 14;
      figure('color','w','name',fldnm{f});
      [ax,h1,h2] = plotyy(1:length(mu_S),mu_S,1:length(mu_S)-1,MAd_mu_MA);
      set(h1,'linestyle','none','marker','.','displayname','Fixed Horizontal Avg ($\mu_L%)');
      set(h2,'linewidth',2,'linestyle','-','color','c','marker','s',...
        'markersize',3,'markerfacecolor','c',...
        'displayname','Filtered $\frac{d\mu_{MA}}{dx}$ $\left($\frac{d\mu_{MA}}{dx}$\right)_{MA}$');
      hold(ax(1), 'on');
      ylabel(ax(1),'Intensity (cps)','fontsize',fntsz);
      plot(ax(1),1:length(mu_MA),mu_MA,'.g',...
        'displayname','Filtered $\mu_S$ ($\mu_{MA}$)');
      %plot(ax(1),1:length(mu_MA),mean(v_non0,2)','*',...
      %  'displayname','Original Mean','color',[100,0,200]./256);
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
      set(ax(2),'xgrid','on','ygrid','on');
      xlabel(ax(2),'Line','fontsize',fntsz);
      ylabel(ax(2),'Intensity (cps)','fontsize',fntsz);
      set(ax,'xlim',[0,length(mu_S)])
      yl1= get(ax(1),'ylim');
      yl2 = get(ax(2),'ylim');
      plot(ax(2),[up_co,up_co],yl2,'--r','linewidth',2);
      text(up_co+2,yl1(2)-(yl1(2)-yl1(1))/5,num2str(up_co),'fontsize',fntsz);
      plot(ax(2),[low_co,low_co],yl2,'--r','linewidth',2);
      text(low_co+2,yl1(2)-(yl1(2)-yl1(1))/5,num2str(low_co),'fontsize',fntsz);
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