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
%  23 Sept 2015  Amanda G Balderrama   amandagbalderrama@gmail.com    5
%     Utilized newer functions for region labeling
%  23 Sept 2015  Amanda G Balderrama   amandagbalderrama@gmail.com    5.1
%     Rather than smoothing a line of data, data are spatially smoothed and
%     then statistics are computed
%  17 Feb  2017  Amanda G Balderrama   amandagbalderrama@gmail.com    6
%     Fixed for compatiblity with newest versions of functions


PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

dnorm = [];
stat = 1;

if strmatch('hdrtxt',PropertyNames)
  hdrtxt = PropertyVal{strmatch('hdrtxt',PropertyNames)};
end

if strmatch('typestr',PropertyNames)
  typestr = PropertyVal{strmatch('typestr',PropertyNames)};
else
  typestr = '';
end

if strmatch('labelinfo',PropertyNames)
  labelInfo = PropertyVal{strmatch('labelinfo',PropertyNames)};
else
  labelInfo = '';
end

[data,hdrtxt,d_div,line_types,type_rows] = divMIMSstruct(data,hdrtxt);
if isempty(labelInfo)
  [data,hdrtxt,d_div,line_types,type_rows,labelInfo] = labelMIMSstruct(data,hdrtxt);
end

%[data,hdrtxt,d_div,line_types,type_rows,labelMask,laONmask,tsumask,dchMask] = ...
%  labelMIMSstruct(data,hdrtxt,'ncc',1);

if isempty(typestr)
  typestr = line_types;
elseif isstr(typestr)
  typestr = {typestr};
end

if isempty(elem_ana_range)
  elem_ana_range = fieldnames(rmfield(data,skip_fields));
end
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
nel = length(elem_ana_range);

for c = 1:length(typestr)
  typec = typestr{c};
  usetype = find(cellfun(@(x) ~isempty(x),strfind(line_types,typec)));
  dtype = d_div{usetype}; %dtype.dataset = evalf{1};
  typechMask = labelInfo.dchMask(type_rows{usetype},1:size(labelInfo.labelblocks{usetype},2),:);
  if any(any(typechMask(:,:,4)))
    tsumask = imdilate(~isnan(typechMask(:,:,4)),strel('octagon',6));
    tmaskprop = regionprops(tsumask);
    bkgnd_ind = [1:floor(tmaskprop.BoundingBox(1))];%,... To incorporate both before and after data, uncomment the following line. 
      %ceil(tmaskprop.BoundingBox(1)+tmaskprop.BoundingBox(3)):size(dtype.Time,2)];
    la0ind = find(sum(sum(~isnan(typechMask(:,:,1:2)),3)) > size(dtype.Time,2)*0.01);
    usebki = find(~ismember(bkgnd_ind,la0ind));
    bkgnd_ind = bkgnd_ind(usebki);
  end
  
  for el = 1:nel
    if isempty(strmatch(fldnm{el},'Time')) && isempty(strmatch(fldnm{el},'line_datevec')) && el <= length(fldnm)
      dtel = getfield(dtype,fldnm{el});
      [M,N] = size(dtel);
      dnan = isnan(dtel);
      dnani = find(dnan);
      for i = 1:length(dnani)
        [ri,ci] = ind2sub(size(dtel),dnani(i));
        dtel(ri,ci) = nanmean(dtel(ri,max([ci-3,1]):min([ci+3,N])));
      end
      included_rows = 1:M;
      vbkgd = dtel(:,bkgnd_ind);
      
      %       [~,bN] = size(vbkgd);
      %       filth = 3;
      %       filtw = round(bN/10);
      %       %vbkgd = medfilt2(vbkgd,[filth,filtw],'symmetric');
      %       uval = unique(vbkgd);
      %       luval = log(uval);
      %       [f,xi,bw] = ksdensity(log(vbkgd(vbkgd ~=0)),'bandwidth',0.1);
      %       hxi = hist(log(vbkgd(vbkgd ~=0)),xi);
      %       figure; plot(xi,f);hold on;
      %       plot(xi,hxi./max(hxi).*max(f),'c');
      %       figure; imagesc(log(vbkgd(:,30:end)),[xi(1),xi(end)]);
      %       figure; hold on;
      %       linec = colormap(jet(M));
      %       for i = 3:30:M-3
      %         rdat = log(vbkgd((i-2):(i+2),30:end));
      %         uh = ksdensity(rdat(:),xi);
      %         plot(xi(uh~=0),uh(uh~=0),'color',linec(i,:),'displayname',sprintf('row %d-%d',i-2,i+2));
      %       end
      
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
        legend({'FT\{mean\}','FT\{median\}'},'fontsize',fntsz);
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
        switch norm_dat
          case 'mean'
            norm_dat = line_mean;
            norm_str = 'Mean';
          case 'med'
            norm_dat = line_med;
            norm_str = 'Median';
          case 'lpfmean'
            norm_dat = lpfmean;
            norm_str = 'LPF Mean';
          case 'lpfmed'
            norm_dat = lpfmed;
            norm_str = 'LPF Median';
        end
      else
        norm_dat = lpfmed;
        norm_str = 'LPF Median';
      end
      
      dnon0_norm = dtel./repmat(norm_dat,1,N);  %<== commented
      %out 29 Aug because of testing on 120710Brain 'Mo95' data. Multiple
      %datum had intensity values of 0 resulting in huge numbers in the
      %division. Method instead subtracts mean from all numbers in a given
      %line, then ensures there are no 0 values.
      %
      % QUESTION REMAINS AS TO WHETHER line_mean, line_median, LPFlnmed, or
      % LPFlnmean should be subtracted.
      %dnon0_norm = d - repmat(norm_dat,1,size(d,2));
      if min(dnon0_norm(:)) < 0
        dnon0_norm = dnon0_norm + abs(min(dnon0_norm(:)));
      end
      d_norm = nan(size(dtel));
      d_norm(included_rows,:) = dnon0_norm;
      if strmatch('heat',PropertyNames)
        heatmap = PropertyVal{strmatch('heat',PropertyNames)};
        if heatmap
          figure('position',[94         445        1245         408]);
          dth = auto_thresh(dtel);
          h1 = subplot(1,2,1); 
          plot_heatmap(dtel,sprintf('%s',fldnm{el}),...
            'thresh',dth,'plot','lin','handle',h1);
          h2 = subplot(1,2,2);
          plot_heatmap(d_norm,sprintf('%s Normalized by %s',fldnm{el},norm_str),...
            'thresh',dth,'plot','lin','handle',h2);
        end
      end
      dnorm = setfield(dnorm,fldnm{el},d_norm);
      
      if strmatch('stat',PropertyNames)
        stat = PropertyVal{strmatch('stat',PropertyNames)};
      end
      
      if stat
        fntsz = 14;
        figure('color','w'); plot(included_rows,[line_mean,line_med,lpfmean,lpfmed,...
          mean(d_norm(included_rows,bkgnd_ind),2),...
          median(d_norm(included_rows,bkgnd_ind),2)]);%,'-*');
        set(gca,'ygrid','on','fontsize',fntsz);
        alld = [line_mean,line_med,lpfmean,lpfmed];
        %ylim([min(min(alld)), max(median(d_norm(included_rows,bkgnd_ind),2))*1.1])
        xlabel('Line Number','fontsize',fntsz);
        ylabel('Intensity (cps)','fontsize',fntsz);
        legend('mean','median','LPFmean','LPFmedian','d_{norm} mean','d_{norm} median',...
          'location','southeast');
        title(sprintf('Background metrics for %s -- Normalized by %s',fldnm{el},norm_str),...
          'fontsize',fntsz);
      end
      
    elseif strmatch(fldnm{el},'Time')
      dnorm = setfield(dnorm,'Time',dtel.Time);
    elseif strmatch(fldnm{el},'line_datevec')
      dnorm = setfield(dnorm,'line_datevec',dtel.line_datevec);
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
  
  vco_out = [];
end