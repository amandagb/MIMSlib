function [phi_MI,varargout] = twoimgMI(A,B,varargin)
%% Script: [phi_MI] = twoimgMI(A,B,varargin)
% Description: The following implementation is derived from 
%   [Mattes et al.(2003)] D. Mattes, D. Haynor, H. Vesselle, T. Lewellen,
%   and W. Eubank. PET-CT image registration in the chest using free-form
%   deformations. Medical Imaging, IEEE Transactions on, 22(1):120 –128,
%   Jan 2003. ISSN 0278-0062. doi: 10.1109/ TMI.2003.809072.
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% A:  Image 1 (fixed/reference image)
% B:  Image 2 (moving/test image)
% varargin - 'PropertyName','PropertyValue'
%   'mu':    1x6 vector with transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations [DEFAULT = 50]
%   'klog_cell': 1 x nbins cell with each cell containing the M x N (image
%   dimensions of A) 0-spline values from the fixed image, A, histogram
%   'pA': 1 x nbins vector which approximates the pdf of the fixed image, A
%   'plot': binary indicating whether to plot (1) or not (0) [DEFAULT = 0]
%   'axh': axis handle which should be used for plotting
% OUTPUTS ---------------------------------------------------------------
% phi_MI: Mutual information between fixed image A and moving image B
% varargout
%   {1} structure with fields which contain the distributions for A (pA), B
%   (pB) and the joint (pAB)
% 
%
%  Date           Author            E-mail                      Version
%  22 Oct 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  23 Oct 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  09 Apr 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%  02 May 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Found that after transforming B (that is variable TxB), the values in
%     TxB could be outside of the range of the values in B. This was fixed
%     by adding lines 86,87
%  22 July 2013  Amanda Gaudreau   amanda.gaudreau@gmail.com     4.1
%     Added axis numbering such that it reflects the original image
%     intensities rather than the bin numbers
%  6  Sept 2013  Amanda Gaudreau   amanda.gaudreau@gmail.com     4.2
%     Adding data structures to mirror time points recorded in twoimgMImd

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
  if p
    if strmatch('axh',PropertyNames)
      axh = PropertyVal{strmatch('axh',PropertyNames)};
    else
      figure;
      axh = get(gcf,'currentaxes');
    end
  end
else
  p = 0;
end

if strmatch('mu',PropertyNames)
  mu = PropertyVal{strmatch('mu',PropertyNames)};
else
  mu = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk
end

if strmatch('logbase',PropertyNames)
  logbase = PropertyVal{strmatch('logbase',PropertyNames)};
else
  logbase = '2';
end

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 50;
end
[Mr,Nc] = size(A);
Abins = nbins;
Bbins = nbins;

% ------------------------------------------------------------------------
% • Implements the zero order B-spline definition in the for loop to
% represent the vector B image as bin indices rather than original
% intensity values -- these will be used to guide the relevant bins for the
% 3rd order B-spline
% ------------------------------------------------------------------------
w = 1; c(w,:) = clock; % w = 1
w = w + 1;

if strmatch('klog_cell',PropertyNames)
  klog_cell = PropertyVal{strmatch('klog_cell',PropertyNames)};
else
  klog_cell = [];
  
  Arng = max(A(:)) - min(A(:));
  bin_rngA = Arng/(Abins-1); % equally spaced bins which cover range of A
  Aind = (A-min(A(:)))./bin_rngA; % index image which reference which bin the pixel belongs to
  Alims = [0.5:1:(Abins-0.5)].*bin_rngA + min(A(:));
end

if strmatch('pA',PropertyNames)
  pA = PropertyVal{strmatch('pA',PropertyNames)};
else
  pA = [];
end

if isempty(klog_cell);
  klog_cell = cell(Abins,1);
  [klog_cell{1:end}] = deal(zeros(Mr,Nc));
  for k = 1:Abins
    b0arg = (k-1) - Aind;
    klog_cell{k} = (b0arg > -1/2 & b0arg <= 1/2); %find(b0arg <= 0 & b0arg > -1);
  end
end
% 
% if isempty(pA); 
%   pA = zeros(1,Abins); 
%   for k = 1:Abins
%     pA(k) = sum(klog_cell{k}(:));
%   end
%   normfac = (sum(pA));
%   pA = pA./normfac;
% end

c(w,:) = clock; % w = 2
Aformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • 
% ------------------------------------------------------------------------
TxB = transform_image(B,mu,'extrapval',nan,'outputimsz',size(A));%,'interpmethod','nearest'); % added property 'outputimsz' 2013 April 9
nanmat = isnan(TxB); %disp(sum(nanmat))
nonnanID = nanmat == 0;%find(nanmat == 0);
TxB(TxB > max(B(:))) = max(B(:));
TxB(TxB < min(B(:))) = min(B(:));
Brng = max(B(:)) - min(B(:));
Bbins = nbins;
bin_rngB = Brng/(Bbins-1);
Bind = (TxB-min(TxB(:)))./bin_rngB;
Blims = [0.5:1:(Bbins-0.5)].*bin_rngB + min(B(:));

b3fun_cell = cell(Bbins,1);
[b3fun_cell{1:end}] = deal(zeros(Mr,Nc));
for l = 1:Bbins
  b3arg = abs((l-1) - Bind);
  l0to1logical = (b3arg <= 1);
  l1to2logical = (b3arg < 2 & b3arg > 1);
  b3fun_cell{l}(l0to1logical) = 1/6.*(4 - 6.*b3arg(l0to1logical).^2 + 3.*b3arg(l0to1logical).^3);
  b3fun_cell{l}(l1to2logical) = 1/6.*(2-b3arg(l1to2logical)).^3;
end

c(w,:) = clock; % w = 3
Bformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

c(w,:) = clock; % w = 4
unqABpairtime = nan;%etime(c(w,:),c(w-1,:));
w = w + 1;
c(w,:) = clock; % w = 5
totABpairtime = nan;%etime(c(w,:),c(w-2,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Loop which forms pAB
% ------------------------------------------------------------------------
pAB = zeros(Abins,Bbins);
for k = 1:Abins
  for l = 1:Bbins
    non0k = find(klog_cell{k} ~= 0);
    usebind = nonnanID(non0k) == 1;
    kl_img = b3fun_cell{l}(non0k(usebind));
    pAB(k,l) = sum(kl_img); % sum(Bind(non0k(usebind)) >= (l - 1/2) & ...
      %Bind(non0k(usebind)) < (l + 1/2));%
  end
end
c(w,:) = clock; % w = 6
ABpairlooptime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Only perform the following operations on non-zero indices since all
% others will not contribute to the final result
% ------------------------------------------------------------------------
non0pAB = find(pAB ~= 0);

c(w,:) = clock; % w = 6
non0pABtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Normalizing pAB
% ------------------------------------------------------------------------
pAB(non0pAB) = pAB(non0pAB)./sum(pAB(non0pAB));

c(w,:) = clock; % w = 7
pABscaletime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Create pA by summing along dA dimensions
% ------------------------------------------------------------------------
pA = sum(pAB,2);

c(w,:) = clock; % w = 8
pAcalctime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Create pB by summing along dA dimensions
% ------------------------------------------------------------------------
pB = sum(pAB,1);

c(w,:) = clock; % w = 9
pBcalctime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute denominator of MI log function 
% ------------------------------------------------------------------------
pkpl_product =  pA(:)*pB(:)';%repmat(pA(:),1,length(pB)).*repmat(pB(:)',length(pA),1);

c(w,:) = clock; % w = 10
pApBprodtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute MI argument for each relevant element in distriubtions
% ------------------------------------------------------------------------
if strmatch(logbase,'2')
  MIarg = pAB(non0pAB).*log2((pAB(non0pAB))./(pkpl_product(non0pAB)));
else
%   MIarg = pAB.*log((pAB+eps)./(pkpl_product+eps));
  MIarg = pAB(non0pAB).*log((pAB(non0pAB))./(pkpl_product(non0pAB)));
end

c(w,:) = clock; % w = 11
MIargtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Actual computation of MI by summing the MI arguments generated before
% ------------------------------------------------------------------------
phi_MI = sum(MIarg(:));

c(w,:) = clock; % w = 11
MIargsumtime = etime(c(w,:),c(w-1,:));

tot_time = etime(c(w,:),c(1,:));

all_times.Aformtime = Aformtime;
all_times.Bformtime = Bformtime;
all_times.unqABpairtime = unqABpairtime;
all_times.totABpairtime = totABpairtime;
all_times.ABpairlooptime = ABpairlooptime;
all_times.non0pABtime = non0pABtime;
all_times.pABscaletime = pABscaletime;
all_times.pAcalctime = pAcalctime;
all_times.pBcalctime = pBcalctime;
all_times.pApBprodtime = pApBprodtime;
all_times.MIargtime = MIargtime;
all_times.MIargsumtime = MIargsumtime;
all_times.tot_time = tot_time;

dist.pA = pA;
dist.Alims = Alims;
dist.pB = pB;
dist.Blims = Blims;
dist.pAB = pAB;
dist.time = all_times;
varargout{1} = dist;

% f1 = [15         610         500         375];
% f2 = [495         554         500         375];
% f3 = [495         610         500         375];
% f4 = [968         610         500         375];
% f5 = [1440         610         500         375];
% figure('position',f1); imagesc(TxB); colormap(gray(50)); colorbar; title(TxB);
% figure('position',f2); hist(Bind(:),50); axis tight;xlabel('Bind values')
% figure('position',f3); hist(TxB(:),50); axis tight; xlabel('TxB values')
% figure('position',f4); bar(pB); axis tight; title('pB')
% figure('position',f5); imagesc(pAB); colorbar;xlabel('B'); ylabel('A')

if p
  set(gcf,'currentaxes',axh)
  ax1 = subplot('position',[0.2,0.2,0.7,0.7]); imagesc(pAB);
  set(gca,'xticklabel','','yticklabel','')
  subplot('position',[0.1,0.2,0.1,0.7]); plot(pA,1:nbins,'marker','.'); axis ij;
  hold on; %plot(sum(pAB,2),1:nbins,':r');
  set(gca,'xticklabel',''); hold off;
  %-----------------------------------------------------
  % comment out these lines to use bin numbers rather than intensities
  Aticks = str2num(get(gca,'ytickLabel'));
  ntick = length(Aticks);
  ntick = 11;
  tickdist = round((max(A(:)) - min(A(:)))/(ntick-1));
  Aticklab = floor(min(A(:))):tickdist:ceil(max(A(:)));
  Aticklab(1) = ceil(min(A(:))*10)/10;
  Aticklab(end) = floor(max(A(:))*10)/10;
  tickpos = (Aticklab-min(A(:)))./bin_rngA; % index image which reference which bin the pixel belongs to
  set(gca,'ytick',tickpos,'yticklabel',Aticklab);
  %useAticks = Aticks(Aticks > 0);
  %set(gca,'ytick',useAticks,'yticklabel',Alims(useAticks));
  ylabel('Fixed image intensities, $$x$$','interpreter','latex')
  xlabel('$$p(x)$$','interpreter','latex')
  %-----------------------------------------------------
  subplot('position',[0.2,0.1,0.7,0.1]); plot(1:nbins,pB,'marker','.');
  hold on; %plot(1:nbins,sum(pAB,1),':r');
  axis tight;
  set(gca,'yticklabel',''); hold off;
  %-----------------------------------------------------
  % comment out these lines to use bin numbers rather than intensities
  Bticks = str2num(get(gca,'xtickLabel'));
  ntick = length(Bticks);
  ntick = 11;
  tickdist = round((max(TxB(:)) - min(TxB(:)))/(ntick-1));
  Bticklab = floor(min(TxB(:))):tickdist:ceil(max(TxB(:)));
  Bticklab(1) = ceil(min(TxB(:))*10)/10;
  Bticklab(end) = floor(max(TxB(:))*10)/10;
  tickpos = (Bticklab-min(TxB(:)))./bin_rngB; % index image which reference which bin the pixel belongs to
  set(gca,'xtick',tickpos,'xticklabel',Bticklab);
  %useBticks = Bticks(Bticks > 0);
  %set(gca,'xtick',useBticks,'xticklabel',Blims(useBticks));
  xlabel('Moving image intensities, $$y$$','interpreter','latex')
  ylabel('$$p(y)$$','interpreter','latex')
  %-----------------------------------------------------
  %   subplot('position',[0.1,0.1,0.1,0.1]);
  %   imagesc(pABindep);
  %   set(gca,'xticklabel','','yticklabel','')
  set(gcf,'currentaxes',ax1)
end
varargout{1} = dist;
end