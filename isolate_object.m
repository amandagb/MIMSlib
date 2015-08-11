function [mask,opt_thresh,hind,vind,normparams] = isolate_object(I,varargin)
%% Script:  [mask, opt_thresh,hind,vind,normparams] = isolate_object(MIcu_PIdim);
% Description:  Given an input image, finds regions on interest (which are
% determined by modality of the histogram -- typically bimodal). There are
% several approaches which could be used to achieve this. The one
% implemented in v1 is a hypothesis testing threshold on the histogram to
% separate a bimodal distribution into H0 (background) and H1 (foreground).
% Another strategy would be to perform edge detection and isolate
% background and foreground based on edges.
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% I:    Bimodal data with Gaussian kernel
% varargin - 'PropertyName','PropertyValue'
%   'plot': binary indicating whether to display plots (1) or not (0)
%   'nbins': number indicating the number of bins that should be used for
%     the histogram [DEFAULT = 1/15 of max value]
%   'thresh': binary indicating whether to threshold the data (1) or not (0)
%     [DEFAULT = 0]. The property value can also be a 1x2 vector
%   'lpfhist': binary indicating whether to lpf the histogram or not
% OUTPUTS ---------------------------------------------------------------
% mask: binary image which results from thresholding the intensities in I
%   using the opt_thresh
% opt_thresh: the intensity value which corresponds to the optimal division
%   between the two distributions
% hind: horizontal indices which contain the object
% vind: verticle indices which contain the object
% normparams: # modes x 2 matrix where each row corresponds to a
%   distribution and the first column corresponds to the position of the
%   mean and the second column corresonds to the standard deviation
%
%  Date           Author            E-mail                      Version
%  27 June 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  24 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Uses a partial histogram, progressive estimation method to solve for
%     the standard deviation (see estimate_sigma sub-function description)
%  25 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Will implement a MMSE (?) method for fitting data to the sum of two
%     histograms
%  30 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Different method of histogramming

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end

I = double(I);
Ivec = I(:);
if size(I,3) == 3; I = double(rgb2gray(uint8(I))); end
% [Iunique,indvec,bins] = unique(Ivec); %indu gives is the same length as Ivec and each element of the vector corresponds to the index of Iunique
% Ihist = hist(bins,length(Iunique));

if strmatch('thresh',PropertyNames)
  thresh = PropertyVal{strmatch('thresh',PropertyNames)};
else
  thresh = 0;
end

if thresh
  if thresh == 1
    th = auto_thresh(Ivec,'auto',[]);
  else
    th = auto_thresh(Ivec,thresh,[]);
  end
  MNorig = length(Ivec);
  Ivec = Ivec(Ivec <= th(2));
  Ivec = Ivec(Ivec >= th(1));
  MNth = length(Ivec);
end

if strmatch('nbin',PropertyNames)
  nbins = PropertyVal{strmatch('nbin',PropertyNames)};
else
  nbins = round(length(unique(Ivec))/15);
end

[Ihist,xhist] = hist(Ivec,nbins);
if strmatch('lpfhist',PropertyNames)
  lpfhist = PropertyVal{strmatch('lpfhist',PropertyNames)};
else
  lpfhist = 0;
end

if lpfhist
  Ihist_orig = Ihist;
  Ihist = lpf_data(Ihist,16,0.07);
else
  Ihist_orig = Ihist;
end

unqval = unique(Ihist);
sortval = sort(unqval);
rng = sortval(end) - sortval(1);
startv = mean(Ihist)/2;
[v,nearesti] = min(abs(sortval-startv));
bottom_thresh = 0;
i = nearesti;
if p; figure; end
while ~bottom_thresh 
  binaryvec = Ihist > sortval(i);
  changevec = diff(binaryvec);
  ngrps = sum(abs(changevec));
  if p
    if lpfhist; plot(xhist,Ihist_orig,'g'); hold on; end
    plot(xhist,Ihist,'marker','.'); hold on;
    plot([min(Ivec),max(Ivec)],sortval(i).*ones(1,2),'k:');
    plot(xhist,rng*binaryvec,'r');
    plot(xhist(2:end),rng/2*(1+changevec),'m*');
    legend('Original Hist','Binary Thresh','Hist > Thresh','Binary Changes');
    title(sprintf('i = %d, v = %g',i,sortval(i)));
  end
  if ngrps > 1 && rem(ngrps,2) < 2 && binaryvec(1) == 1
    bottom_thresh = sortval(i);
    changevec(1) = 1;
  elseif rem(ngrps,2) == 0 && ngrps > 2
    bottom_thresh = sortval(i);
  else
    i = i + 1;
  end
end

startchng = find(changevec == 1); 
endchng = find(changevec == -1)-1;
nmodes = length(startchng);
bin_ind = ones(nmodes,2);
bin_ind(:,2) = length(Ihist);

for i = 1:nmodes
  irng = startchng(i):endchng(i);
  [mode_max(i),modei] = max(Ihist(irng));
  m(i) = xhist(irng(modei)); %bins(irng(modei));
  if i > 1
    midpnt(i-1) = round((m(i) - m(i-1))/2 + m(i-1));
    [v,bin_ind(i-1,2)] = min(abs(midpnt(i-1) - xhist));
    bin_ind(i,1) = bin_ind(i-1,2) + 1;
    if i == 2
      M{i-1} =  Ihist(bin_ind(i-1,1):bin_ind(i-1,2));%Ivec(Ivec <= Iunique(bin_ind(i-1,2)));% 
      M{i} = Ihist(bin_ind(i,1):bin_ind(i,2));%(Ivec >= Iunique(bin_ind(i,1))); %
    else
      M{i} = Ihist(bin_ind(i,1):bin_ind(i,2));
    end
  end
end
if p
  plot(midpnt(1).*ones(1,2),[min(Ihist),max(Ihist)],'k:'); hold off; 
  ylim([min(Ihist),max(Ihist)])
end
sig_est = zeros(nmodes,2);
fine_sig = -1:0.1:1;
if p; figure; end
for i = 1:length(startchng)
  Mx{i} = xhist(bin_ind(i,1):bin_ind(i,2));
  [sig_est(i,1),sig_est(i,2)] = estimate_sigma(Ivec,M{i},Mx{i},m(i),[1.5:0.1:2.5]);
  if i == 1 && changevec(1) == 1
    s(i) = sig_est(i,2);
  elseif i == 1
    s(i) = sig_est(i,1);
  elseif i == nmodes
    s(i) = sig_est(i,2);
  else
    s(i) = mean([sig_est(i,:)]);
  end
  
%   figure;
%   for n = 1:length(fine_sig)
%     j = fine_sig(n);
%     Ctest = max(M{i})*sqrt(2*pi*(sig(i)+j)^2);
%     Mtstnorm = M{i}./Ctest;
%     pdftest_est = (1/sqrt(2*pi*(sig(i)+j)^2)).*exp(-(Mx{i}-mode_means(i)).^2./(2*(sig(i)+j)^2));
%     fiterr(n) = sum((pdftest_est-Mtstnorm).^2);
%     plot(Mx{i},[Mtstnorm;pdftest_est]); title(sprintf('%d',n));
%   end
%   [v,bestsig] = min(fiterr);
%   sig(i) = sig(i) + fine_sig(bestsig);
  Chist(i) = max(M{i})*sqrt(2*pi*s(i)^2);
  Mnorm{i} = M{i}./Chist(i);
  pdf_est{i} = (1/sqrt(2*pi*s(i)^2)).*exp(-(Mx{i}-m(i)).^2./(2*s(i)^2));
  %Mexp = repmat(0.5:0.1:1.5,length(Mx{i}),1).*repmat(exp(-Mx{i}),1,11).^(repmat(0.5:0.1:1.5,length(Mx{i}),1));
  %Mexp = Mexp./repmat(sum(Mexp),length(Mx{i}),1);
  if p;
    plot(Mx{i},[Mnorm{i}(:),pdf_est{i}(:)],'marker','.');%plot(bin_ind(i,1):bin_ind(i,2),[Mnorm{i}';pdf_est{i}']);% 
    hold on;
  end
  %plot(Mx{i},Mexp','.');
end
normparams = [m',s'];
p0 = 2/3;
p1 = 1 - p0;
eta = p0/p1;
yp = (m(2)*s(1)^2 - m(1)*s(2)^2 + ...
  s(1)*s(2)*sqrt(2*s(2)^2*log((s(2)/s(1))*eta) - 2*s(1)^2*log((s(2)/s(1))*eta) - ...
  2*m(1)*m(2) + m(1)^2 + m(2)^2))/(s(1)^2-s(2)^2);
ym = -(m(1)*s(2)^2 - m(2)*s(1)^2 + ...
  s(1)*s(2)*sqrt(2*s(2)^2*log((s(2)/s(1))*eta) - 2*s(1)^2*log((s(2)/s(1))*eta) - ...
  2*m(1)*m(2) + m(1)^2 + m(2)^2))/(s(1)^2-s(2)^2);
% sig = mean(s);
% g = (m(2)+m(1))/2 + (sig^2*log(eta))/(m(2)-m(1))
opt_thresh = ym;
I10 = zeros(size(I));
I10(I >= opt_thresh) = 1;
SE = strel('rectangle',[6,6]);%('arbitrary',ones(4,4));%[0,1,0;1,1,1;0,1,0]);
mask = imopen(I10,SE);
if p;
  plot_heatmap(I,[],'thresh','auto'); %figure; imshowpair(I10,mask); colormap(gray);
end

hsum = sum(mask,1);
fhsum = lpf_data(hsum,30,0.07);
[thresh_ind,min_changes,interm_thresh,hind] = find_change_thresh(hsum,'starting','b');
hind = thresh_ind(1)+1:thresh_ind(2)-1;
if p;
  hold on; plot(thresh_ind(1).*ones(1,2),[0,size(mask,1)],'y');
  hold on; plot(thresh_ind(2).*ones(1,2),[0,size(mask,1)],'y');
end

vsum = sum(mask,2);
fvsum = lpf_data(vsum,30,0.07);
[thresh_ind,min_changes,interm_thresh,vind] = find_change_thresh(vsum,'starting','b');
vind = thresh_ind(1)+1:thresh_ind(2)-1;
if p;
  hold on; plot([0,size(mask,2)],thresh_ind(1).*ones(1,2),'y');
  hold on; plot([0,size(mask,2)],thresh_ind(2).*ones(1,2),'y');
end

mask = I10;

  function [Lmean,Rmean] = estimate_sigma(data,hist_dat,xbins,mu,sig_weight)
    % Given the histogram follows a Gaussian distribution, it is defined:
    % (1)   H(x)    = C / sqrt(2 pi s^2) * exp(- (x-m)^2 / 2s^2)
    % where C is a constant since H(x) isn't a pdf, m is the mean and s is
    % the standard deviation. H(x) is maximal at the mean of the
    % "distribution", x = m, where 
    % (2)   H(m)    = C / sqrt(2 pi s^2). 
    % Since the distribution is subject to interferences farther from the
    % mean, the for loop looks to the right and to the left of the mean a
    % given number of standard deviations away. For instance, when w = 1 =
    % s (where w is the weight on the standard deviation s), we have:
    %       H(m+s)  = C / sqrt(2 pi s^2) * exp(- (m+s-m)^2 / 2s^2)
    %               = C / sqrt(2 pi s^2) * exp(- 1 / 2) 
    %               = H(m) * exp(- 1 / 2) 
    % When w = 1.5, we have:
    %       H(m+1.5s)= C / sqrt(2 pi s^2) * exp(- (m+1.5s-m)^2 / 2s^2)
    %               = C / sqrt(2 pi s^2) * exp(- (1.5)^2 / 2) 
    %               = H(m) * exp(- 1.5^2 / 2) 
    % So generally, for a given w, we have:
    % (3)   H(m+w*s)= C / sqrt(2 pi s^2) * exp(- (m+w*s-m)^2 / 2s^2)
    %               = H(m) * exp(- w^2 / 2) 
    % Since the data values are found simply using the histogram levels, we
    % can solve for the respective "s" which would give rise to a given
    % histogram value at a certain distance. Knowing the x location where
    % the histogram get to the value of H(m+sig_fram*s) [gtsigvec], we can
    % solve for sigma using the equation:
    % (4)   Lx = m + w*s => s = (Lx - m) / w
    % (5)   Rx = m - w*s => s = (m - Rx) / w
    
    Hm = max(hist_dat); % See (2) in description
    nsig = length(sig_weight);
    Lsig = zeros(1,nsig);
    Rsig = zeros(1,nsig);
    for isig = 1:length(sig_weight)
      HmpWs = Hm*exp(-sig_weight(isig)^2/2); % See (3) in description
      gtsigvec = find(hist_dat >= HmpWs);
      [v,Iind] = min(abs(data - xbins(gtsigvec(1))));
      Lx = data(Iind);
      Lsig(isig) = (mu - Lx)/sig_weight(isig); % See (4) in description
      [v,Iind] = min(abs(data - xbins(gtsigvec(end))));
      Rx = data(Iind);
      Rsig(isig) = (Rx - mu)/sig_weight(isig); % See (5) in description
    end
    Lmean = mean(Lsig);
    Rmean = mean(Rsig);
    avg_sig = mean([Lmean,Rmean]);
  end
end



