function [phi_MI,varargout] = twoimgMIkde(A,B,varargin)
%% Script: [phi_MI] = twoimgMImd(A,B,varargin)
% Description: The following implementation is derived from
%   [Hansen(2009)] B. Hansen. Lecture notes on nonparametrics. ECON 718,
%     2009. URL http://www.ssc.wisc.edu/~bhansen/718/718.htm.
%   Implementation assumes a separable kernel in each dimension. All
%   kernels are second order kernels (nu = 2)
% Example: [MImd,dmd] = twoimgMIkde(F,M,'extrapval',nan,'nbins',binsvec(n),...
%             'plot',0,'logbase','e','Abinlim',Fbinlim,'Bbinlim',Mbinlim,...
%             'kernel',{kernstr{ks}});
% Required Functions: transform_image
% INPUTS ----------------------------------------------------------------
% A:  MrA x NcA Image with dA channels (fixed/reference image). May also be
% a MrA·NcA x dA vector. if NcA is < 10, it will automatically be treated
% as the latter.
% B:  MrB x NcB Image with dB channels (moving/test image) -- see A
% varargin - 'PropertyName','PropertyValue'
%   • 'mu':    1x6 vector with transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   • 'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations [DEFAULT = 50]
%   • 'logbase': string indicating whether the log base should be '2' or 'e'
%   [DEFUALT = 'e']
%   • 'pA': 1 x nbins vector which approximates the pdf of the fixed image, A
%   • 'plot': binary indicating whether to plot (1) or not (0) [DEFAULT = 0]
%   • 'axh': axis handle which should be used for plotting
%   • 'Abinlim': dA x 2 vector indicating the [min,max] bin edge values for
%   each channel
%   • 'Bbinlim': dB x 2 vector indicating the [min,max] bin edge values for
%   each channel
%   • 'kernel': {'string',bandwidth,dimensions} cell identifying and
%   parameterizing the kernel that should be used to construct the pdf
%   estimates. All kernels can be specified or entered as 1-d separable
%   kernels. All kernels are separable second order kernels. The
%   (dA+dB)-dim kernel will be constructed as the product of the separable
%   1-d kernel. Kernels have a functional form K(u) where
%             f_hat(x) = 1/(bw*nsamp)*sum(K((SAMPLES - x)/bw))
%     - 'string': specify the kernel, K(u) => 'norm', 'rect', 'epan', 'b3'
%     [DEFAULT = 'norm']
%     - bandwidth: integer specifying the bandwidth (see f_hat above)
%     [DEFAULT = Silverman's Rule-of-Thumb = sig_hat*C(K)*n^(-1/5)]
%     - dimensions: indicates the dimension the spline should be applied to
%   Multiple kernels can be chosen by letting the variable have multiple
%   rows:
%             {'kernel1',[bw_d1,bw_d2,...],[d1,d2,...];
%             'kernel2',[bw_d1,bw_d2,...],[d1,d2,...]}
%   ..................................................................
%   transform_image varargin
%   ..................................................................
%   • 'extrapval': [DEFAULT = nan]
%   • 'interpmethod': [DEFAULT = 'bicubic']
%   • 'transtype': [DEFAULT = 'physical']
%   • 'fixedstr'  • 'fixedvals'   • 'varstr'    • 'varvals'
%   • 'rotpnt': [DEFAULT = 'center']
%   • 'shape': [DEFAULT = 'same']
%
% OUTPUTS ---------------------------------------------------------------
% phi_MI: Mutual information between fixed image A and moving image B
% varargout
%   {1} structure with fields which contain the distributions for A (pA), B
%   (pB) and the joint (pAB)
%       • pA: distribution of A
%       • Alims: bin centers where values of A range from min(A) to max(A)
%       and bin width can be computed as (Alims(1)-min(A))*2
%       • pB: distribution of B
%       • Blims: bin centers where values of A range from min(A) to max(A)
%       and bin width can be computed as (Alims(1)-min(A))*2
%       • pAB: joint distribution
%
%
%  Date           Author            E-mail                      Version
%  20 Sept 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%     Based on twoimgMImd v3.4
%  27 Sept 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     0.1
%     Attempting to optimize computation to reduce time. This version loops
%     through all possible bins but only evaluates the kernels for
%     non-empty bins.
%  29 Sept 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     0.2
%     Optimize computation to reduce time. Loop only through the closest
%     bins to a given sample and extends the look to bins which fall
%     withing the kernel support (Ksupp). If the support is infinite (ie
%     Gaussian Kernel) then all possible bins in that dimension are
%     evaluated.
%  14 Nov  2016   Amanda Balderrama   amandagbalderrama@gmail.com  1
%     Allowed for abaility to specify the bins in image A or B

[MrA,NcA,dA] = size(A);
allNaNimg = sum(isnan(A),3) < dA;
if NcA == 1 || MrA == 1
  useCols = 1:NcA;
  useRows = 1:MrA;
else
  useCols = find(sum(allNaNimg) > 0);
  useRows = find(sum(allNaNimg') > 0);
end
A = A(useRows,useCols,:);
B = B(useRows,useCols,:);

[MrA,NcA,dA] = size(A);
if NcA < 10
  Aflat = A;
else
  Aflat = reshape(A(:),MrA*NcA,dA);
end
dA = size(Aflat,2);

[MrB,NcB,dB] = size(B);
if NcB < 10
  Bflat = B;
else
  Bflat = reshape(B(:),MrB*NcB,dB);
end
dB = size(Bflat,2);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plot',PropertyNames)
  plot10 = PropertyVal{strmatch('plot',PropertyNames)};
  if plot10
    if strmatch('axh',PropertyNames)
      axh = PropertyVal{strmatch('axh',PropertyNames)};
    else
      figure;
      axh = get(gcf,'currentaxes');
    end
  end
else
  plot10 = 0;
end

if strmatch('mu',PropertyNames)
  mu = PropertyVal{strmatch('mu',PropertyNames)};
else
  mu = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk
end

if strmatch('logbase',PropertyNames)
  logbase = PropertyVal{strmatch('logbase',PropertyNames)};
else
  logbase = 'e';
end

if strmatch('nbinsA',PropertyNames)
  nbinsA = PropertyVal{strmatch('nbinsA',PropertyNames)};
end

if strmatch('nbinsB',PropertyNames)
  nbinsB = PropertyVal{strmatch('nbinsB',PropertyNames)};
end

if strmatch('nbins',PropertyNames)
  if ~exist('nbinsA')
    nbinsA = PropertyVal{strmatch('nbins',PropertyNames)}.*ones(1,dA);
  end
  if ~exist('nbinsB')
    nbinsB = PropertyVal{strmatch('nbins',PropertyNames)}.*ones(1,dB);
  end
else
  nbinsA = 50.*ones(1,dA);
  nbinsB = 50.*ones(1,dB);
end

if strmatch('Abins',PropertyNames)
  Abins = PropertyVal{strmatch('Abins',PropertyNames)};
else 
  Abins = [];
end

if strmatch('Bbins',PropertyNames)
  Bbins = PropertyVal{strmatch('Bbins',PropertyNames)};
  for k = 1:dB
    XevalB{k} = Bbins{k};
    XnumB{k} = [1:1:length(Bbins{k})];
  end
  Bbins = cellfun(@(x) length(x),XevalB);
  if any(Bbins) == 0
    Bbins = [];
    clear XevalB XnumB
  else
    bin_rngB = cellfun(@(x) mean(diff(x)),XevalB);
  end
else 
  Bbins = [];
end

if isempty(Abins)
  if length(nbinsA) < dA
    Abins = nbinsA.*ones(1,dA);
  else
    Abins = nbinsA;
  end
end

if isempty(Bbins)
  if length(nbinsB) < dB
    Bbins = nbinsB.*ones(1,dB);
  else
    Bbins = nbinsB;
  end
end

if strmatch('Abinlim',PropertyNames)
  Abinlim = PropertyVal{strmatch('Abinlim',PropertyNames)};
else
  Abinlim = [];
end

if strmatch('Bbinlim',PropertyNames)
  Bbinlim = PropertyVal{strmatch('Bbinlim',PropertyNames)};
else
  Bbinlim = [];
end

if strmatch('kernel',PropertyNames)
  Kinfo = PropertyVal{strmatch('kernel',PropertyNames)};
else
  Kinfo = {'norm',[]};
end

Kstr = Kinfo(:,1);
if size(Kinfo,2) >= 2
  Kbw = Kinfo(:,2);
else
  Kbw{1} = [];
end

if size(Kinfo,2) >= 3
  Kdim = Kinfo(:,3);
else
  Kdim{1} = 1:dA+dB;
end

% ------------------------------------------------------------------------
% • Turns image A from an [m x n x d] image into an [mn x d] vector
% • Defines the bin centers
% ------------------------------------------------------------------------
w = 1; c(w,:) = clock; % w = 1
w = w + 1;

if isempty(Abinlim)
  minA = min(Aflat);
  maxA = max(Aflat);
else
  minA = Abinlim(:,1)';
  maxA = Abinlim(:,2)';
end
Arng = maxA - minA;
if ~exist('XnumA','var')
  bin_rngA = Arng./Abins;
  XevalA = cell(1,dA); XnumA = cell(1,dA);
  for k = 1:dA
    Aflat(Aflat(:,k) < minA(k),k) = nan;
    Aflat(Aflat(:,k) > maxA(k),k) = nan;
    XevalA{k} = [0.5:1:(Abins(k) - 0.5)].*bin_rngA(k) + minA(k);
    XnumA{k} = [1:1:Abins(k)];
  end
end

c(w,:) = clock; % w = 2
Aformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Turns image B from an [m x n x d] image into an [mn x d] vector
% • Transforms B according to indicate transformation parameters mu
% • Defines the bin centers
% ------------------------------------------------------------------------
TxB = transform_image(B,mu,'extrapval',nan);%,'interpmethod','nearest'); % added property 'outputimsz' 2013 April 9
if NcB < 10
  TxBflat = TxB;
else
  TxBflat = reshape(TxB(:),MrB*NcB,dB);
end
if isempty(Bbinlim)
  maxB = max(Bflat);
  minB = min(Bflat);
else
  minB = Bbinlim(:,1)';
  maxB = Bbinlim(:,2)';
end
Brng = maxB - minB;
Blims = cell(1,dB);
if ~exist('XnumB','var')
  bin_rngB = Brng./Bbins;
  XevalB = cell(1,dB); XnumB = cell(1,dB);
  for k = 1:dB
    TxBflat(TxBflat(:,k) > maxB(k),k) = nan;%maxB(k);
    TxBflat(TxBflat(:,k) < minB(k),k) = nan;%minB(k);
    XevalB{k} = [0.5:1:(Bbins(k)-0.5)].*bin_rngB(k) + minB(k);
    XnumB{k} = [1:1:Bbins(k)];
  end
end

c(w,:) = clock; % w = 3
Bformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Defines joint sample points Xi
% • Defines function evaluation points Xeval
% ------------------------------------------------------------------------
Xi = [Aflat, TxBflat];
useXi = find(sum(isnan(Xi),2) == 0);
Xi = Xi(useXi,:);
npnts = length(useXi);
Xevalcell = [XevalA,XevalB];
Xnumind = [XnumA,XnumB];

% ------------------------------------------------------------------------
% • Parameterize the KDE kernel K -- most parameters identified from Hansen
% 2009 "Lecture Notes on Nonparametrics"
%   nu: kernel order
%   RK: kernel roughness
%   kap2: second moment of the kernel
%     - Parameters are given on table 1, page 4 of the cited document
% ------------------------------------------------------------------------
nu = 2;
Kcell = cell(dA+dB,1);
bw = zeros(1,dA+dB);
Ksupp = nan(2,dA+dB); % support of the kernel
for kd = 1:size(Kdim,1)
  switch lower(Kstr{kd})
    case 'norm'
      RK = 1/(2*sqrt(pi));
      kap2 = 1;
      K = @(u) exp(-u.^2./2)./sqrt(2*pi);
      %Ksupp(:,Kdim{kd}) = repmat([-3;3],1,length(Kdim{kd}));
%     case 'normsupp'
%       RK = 1/(2*sqrt(pi));
%       kap2 = 1;
%       K = @(u) exp(-u.^2./2)./sqrt(2*pi);
%       Ksupp(:,Kdim{kd}) = repmat([-1;1],1,length(Kdim{kd}));
    case 'rect'
      RK = 1/2;
      kap2 = 1/3;
      K = @(u) 0.5.*(abs(u) <= 1);
      Ksupp(:,Kdim{kd}) = repmat([-1;1],1,length(Kdim{kd}));
    case 'epan'
      RK = 3/5;
      kap2 = 1/5;
      K = @(u) 0.75.*(1-u.^2).*(abs(u) <= 1);
      Ksupp(:,Kdim{kd}) = repmat([-1;1],1,length(Kdim{kd}));
    case 'b3'
      RK = 1/3;
      kap2 = 151/315;
      K = @(u) (1/6).*((4-6*u.^2+3.*abs(u)).*(abs(u) < 1) + (2-abs(u)).^3.*(abs(u) >=1 & abs(u) < 2));
      Ksupp(:,Kdim{kd}) =  repmat([-2;2],1,length(Kdim{kd}));
  end
  
  for dim = 1:length(Kdim{kd})
    Kcell{Kdim{kd}(dim)} = K;
  end
  
  if isempty(Kbw{kd})
    sig_hat = std(Xi);
    CK = 2*((pi^(1/2)*factorial(nu)^3*RK)/(2*nu*factorial(2*nu)*kap2^2))^(1/(2*nu+1));
    bw = sig_hat.*CK.*npnts^(-1/(2*nu+1));
  else
    bwkd = Kbw{kd};
    if length(bwkd) < length(Kdim{kd})
      bw(Kdim{kd}) = repmat(bwkd,1,length(Kdim{kd}));
    else
      bw(Kdim{kd}) = bwkd;
    end
  end
end

Ksupp = Ksupp.*repmat(bw,2,1);
c(w,:) = clock; % w = 5
Kformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute pAB using KDE and looping through each bin = 1/(nh) sum_{i=1}^n K((X_i - x)/h)
%  *** NOTE that because pAB is normalized to equal 1 in the end, the extra
%  normalization steps (ie dividing by the kernel BW and # points) is not
%  necessary and only results in an increased computation time. Therefore,
%  these two operations were commented out
% ------------------------------------------------------------------------
binsvec = [Abins,Bbins];
bin_rngs = [bin_rngA,bin_rngB];
pAB = zeros(binsvec); %zeros(prod(binsvec),1);

if sum(isnan(Ksupp(2,:))) == dA+dB
  j = 1;
  ABj = combvec(Xnumind{j});
  Xeval = combvec(Xevalcell{j});
  while j < dA+dB
    j = j + 1;
    ABj = combvec(ABj,Xnumind{j});
    Xeval = combvec(Xeval,Xevalcell{j});
  end
  Xeval = Xeval';
  ABj = ABj';
else
  Xi_int = round((Xi - repmat([minA,minB],npnts,1))./repmat(bin_rngs,npnts,1));
  ABjoints = unique(Xi_int,'rows');
  ABjoints(ABjoints == 0) = 1;
  lasti = 0;
  ABj = nan(dA+dB,min([round(size(ABjoints,1)*max(bw)*5),size(ABjoints,1)]));%round(prod(binsvec)*max(bw)*1));
  for i = 1:size(ABjoints,1)
    j = 0; ieval = [];
    while j < dA+dB
      j = j + 1;
      if isnan(Ksupp(2,j))
        if isempty(ieval)
          ieval = combvec(1:binsvec(j));
        else
          ieval = combvec(ieval,1:binsvec(j));
        end
      else
        if isempty(ieval)
          ieval = combvec(ABjoints(i,j) + [-ceil(Ksupp(2,j)/bin_rngs(j)):ceil(Ksupp(2,j)/bin_rngs(j))]);
          ieval = ieval(:,ieval(j,:) > 0 & ieval(j,:) <= binsvec(j));
        else
          ieval = combvec(ieval,ABjoints(i,j) + [-ceil(Ksupp(2,j)/bin_rngs(j)):ceil(Ksupp(2,j)/bin_rngs(j))]);
          ieval = ieval(:,ieval(j,:) > 0 & ieval(j,:) <= binsvec(j));
        end
      end
    end
    ABj(:,lasti+[1:size(ieval,2)]) = ieval;
    lasti = lasti + size(ieval,2);
  end
  
  ABj = ABj(:,sum(isnan(ABj),1) == 0);
  ABj = unique(ABj','rows');
  Xeval = zeros(size(ABj));
  j = 0;
  while j < dA+dB
    j = j + 1;
    Xeval(:,j) = Xevalcell{j}(ABj(:,j));
  end
end

c(w,:) = clock; % w = 4
Xevaltime = etime(c(w,:),c(w-1,:));
w = w + 1;


% for i = 1:size(Xeval,1)
%   dXi = Xi - repmat(Xeval(i,:),npnts,1);
%   ismembK = find(prod(dXi >= repmat(Ksupp(1,:),npnts,1),2) .* ...
%     prod(dXi <= repmat(Ksupp(2,:),npnts,1),2));
%   %   ismembK = 1:npnts;
%   if ~isempty(ismembK)
%     Kout = zeros(length(ismembK),dA+dB);
%     npnteval = length(ismembK);
%     for kd = 1:size(Kinfo,1)
%       datad = Kdim{kd};
%       Kout(:,datad) = (Kcell{datad(1)}(dXi(ismembK,datad)...
%         ./repmat(h(datad),npnteval,1)));%./repmat(h(datad),npnts,1);
%     end
%     pAB(i) = sum(  prod( Kout,2 ) ) ;%/ (npnts);
%   else
%     pAB(i) = 0;
%   end
% end

for i = 1:size(ABj,1)
  linind = custom_sub2ind(size(pAB),ABj(i,:)); %sum((ABj(i,:)-1).*binsvec.^(0:(dA+dB-1)),2)+1;
  dXi = Xi - repmat(Xeval(i,:),npnts,1);
  if any(isnan(Ksupp))
    ismembK = 1:npnts;
  else
    ismembK = find(prod(dXi >= repmat(Ksupp(1,:),npnts,1),2) .* ...
      prod(dXi <= repmat(Ksupp(2,:),npnts,1),2));
  end
  %   ismembK = 1:npnts;
  if ~isempty(ismembK)
    Kout = zeros(length(ismembK),dA+dB);
    npnteval = length(ismembK);
    for kd = 1:size(Kinfo,1)
      datad = Kdim{kd};
      Kout(:,datad) = (Kcell{datad(1)}(dXi(ismembK,datad)...
        ./repmat(bw(datad),npnteval,1)));%./repmat(h(datad),npnts,1);
    end
    pAB(linind) = sum(  prod( Kout,2 ) ) ;%/ (npnts);
  else
    pAB(linind) = 0;
  end
end

c(w,:) = clock; % w = 6
ABpairlooptime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Only perform the following operations on non-zero indices since all
% others will not contribute to the final result
% ------------------------------------------------------------------------
pAB(pAB < eps) = 0;
non0pAB = find(pAB ~= 0);

c(w,:) = clock; % w = 7
non0pABtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Normalizing pAB
% ------------------------------------------------------------------------
pAB(non0pAB) = pAB(non0pAB)./sum(pAB(non0pAB));

c(w,:) = clock; % w = 8
pABscaletime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Create pA by summing along dA dimensions
% ------------------------------------------------------------------------
Asumdim = (dA+dB):-1:(dA+1);
pA = pAB;
for i = Asumdim
  pA = sum(pA,i);
end

c(w,:) = clock; % w = 9
pAcalctime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Create pB by summing along dA dimensions
% ------------------------------------------------------------------------
Bsumdim = 1:1:dA;
pB = pAB;
for i = Bsumdim
  pB = sum(pB,i);
end
pB = reshape(pB,[Bbins,1]);

c(w,:) = clock; % w = 10
pBcalctime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute denominator of MI log function
% ------------------------------------------------------------------------
pkpl_product = pA(:)*pB(:)';%reshape(pA(:)*pB(:)',binsvec);
% Another opportunity to speed up computation is to look only at the pA, pB
% product pairs which will be used in MIarg. That is, only compute
% length(non0pAB) products which correspond to the nonzero indices in pAB

c(w,:) = clock; % w = 11
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

c(w,:) = clock; % w = 12
MIargtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Actual computation of MI by summing the MI arguments generated before
% ------------------------------------------------------------------------
phi_MI = sum(MIarg(:));

c(w,:) = clock; % w = 13
MIargsumtime = etime(c(w,:),c(w-1,:));

tot_time = etime(c(w,:),c(1,:));

all_times.Aformtime = Aformtime;
all_times.Bformtime = Bformtime;
all_times.Xevaltime = Xevaltime;
all_times.Kformtime = Kformtime;
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
dist.Alims = XevalA;
dist.pB = pB;
dist.Blims = XevalB;
dist.pAB = pAB;
dist.time = all_times;
dist.bw = bw;
dist.kerninfo.bw = bw;
dist.kerninfo.str = Kstr;

dist.pABnon0 = pAB(non0pAB);
dist.logarg = (pAB(non0pAB))./(pkpl_product(non0pAB));
dist.nonNANpxls = npnts;
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

if plot10
  set(gcf,'currentaxes',axh);
  if dA + dB == 2
    ax1 = subplot('position',[0.2,0.2,0.7,0.7]); imagesc(pAB); colormap(jet(32))
    set(gca,'xticklabel','','yticklabel','')
    subplot('position',[0.1,0.2,0.1,0.7]); plot(pA,1:Abins,'marker','.'); axis ij;
    hold on; %plot(sum(pAB,2),1:nbins,':r');
    set(gca,'xticklabel',''); hold off;
    %-----------------------------------------------------
    % comment out these lines to use bin numbers rather than intensities
    Aticks = str2num(get(gca,'ytickLabel'));
    Aticks(1) = 1;
    ntick = length(Aticks);
    %ntick = 6;
    %tickdist = round((max(A(:)) - min(A(:)))/(ntick-1));
    Aticklab = XevalA{1}(Aticks);%floor(min(A(:))):tickdist:ceil(max(A(:)));
    %Aticklab(1) = ceil(min(A(:))*10)/10;
    %Aticklab(end) = floor(max(A(:))*10)/10;
    tickpos = Aticks;%(Aticklab-min(A(:)))./bin_rngA; % index image which reference which bin the pixel belongs to
    set(gca,'ytick',tickpos,'yticklabel',Aticklab);
    %useAticks = Aticks(Aticks > 0);
    %set(gca,'ytick',useAticks,'yticklabel',Alims(useAticks));
    ylabel('Fixed image intensities, $$x$$','interpreter','latex')
    xlabel('$$p(x)$$','interpreter','latex')
    %-----------------------------------------------------
    subplot('position',[0.2,0.1,0.7,0.1]); plot(1:Bbins,pB,'marker','.');
    hold on; %plot(1:nbins,sum(pAB,1),':r');
    axis tight;
    set(gca,'yticklabel',''); hold off;
    %-----------------------------------------------------
    % comment out these lines to use bin numbers rather than intensities
    Bticks = str2num(get(gca,'xtickLabel'));
    if Bticks(1) == 0
      Bticks(1) = 1;
    end
    ntick = length(Bticks);
    %ntick = 6;
    %tickdist = round((max(TxB(:)) - min(TxB(:)))/(ntick-1));
    %Bticklab = floor(min(TxB(:))):tickdist:ceil(max(TxB(:)));
    %Bticklab(1) = ceil(min(TxB(:))*10)/10;
    %Bticklab(end) = floor(max(TxB(:))*10)/10;
    Bticklab = XevalB{1}(Bticks);
    tickpos = Bticks;%(Bticklab-min(TxB(:)))./bin_rngB; % index image which reference which bin the pixel belongs to
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
  elseif dA + dB == 3
    figure; hold on;
    cc = custom_colormaps('jet',16);
    colind = floor(pAB./max(pAB(:)).*size(cc,1)) + 1;
    colind(find(colind == size(cc,1) + 1)) = size(cc,1);
    ppbthresh = 1; % points per bin threshold
    
    for i = 1:size(cc,1)
      linind = find(colind == i);
      non0pAB = pAB(linind) ~= 0 & pAB(linind) > ppbthresh/numel(pAB);
      linind = linind(non0pAB);
      plot3(Xevalcell{1}(pABind(linind,1)),Xevalcell{2}(pABind(linind,2)),Xevalcell{3}(pABind(linind,3)),'.',...'markersize',10,...
        'markerfacecolor',cc(i,:),...
        'markeredgecolor',cc(i,:));
    end
    axis tight;
    grid on;
    xlabel('$$a_1$$','interpreter','latex')
    ylabel('$$a_2$$','interpreter','latex')
    zlabel('$$b_1$$','interpreter','latex')
    %       hold on;
    %       view([40,2])
    %       plot3([mu(G1ind(1)),mu(G1ind(1))],[XevalA{2}(1),mu(G1ind(2))],[XevalB{1}(1),mu(G2ind)],'--r','linewidth',5)
    %       plot3([XevalA{1}(1),mu(G1ind(1))],[mu(G1ind(2)),mu(G1ind(2))],[XevalB{1}(1),mu(G2ind)],'--r','linewidth',5)
    %       plot3([XevalA{1}(1),mu(G1ind(1))],[XevalA{2}(1),mu(G1ind(2))],[mu(G2ind),mu(G2ind)],'--r','linewidth',5)
    %
    %       text(XevalA{1}(2),XevalA{2}(2),XevalB{1}(round(binsvec(bvecind(n))*.95))...
    %         ,{[sprintf('$$\\begin{array}{l} \\mu_1 = [%d, %d] \\\\ \\mu_2 = [%d]\\end{array}$$',mu(G1ind),mu(G2ind))],...
    %         [sprintf('$$\\Sigma = \\left[\\begin{array}{ccc} %1.1f & %1.1f & %1.1f \\\\ %1.1f & %1.1f & %1.1f \\\\%1.1f & %1.1f & %1.1f \\end{array}\\right]$$',Sig_j)]},...
    %         'interpreter','latex','VerticalAlignment','top');
    
    %     title(sprintf('$$p(a_1,a_2,b_1)$$: %d bins, %d $$\\times 10^3$$ samples, $$\\pm  %1.1f \\sigma$$',binsvec(bvecind(n)),tot_samp/1e3,nsig),...
    %       'interpreter','latex');
    
    if dA == 2
      ylabel('$$a_2$$','interpreter','latex')
      zlabel('$$b_1$$','interpreter','latex')
      
      figure; imagesc(XevalA{1},XevalA{2},pA);
      xlabel('$$a_1$$','interpreter','latex')
      ylabel('$$a_2$$','interpreter','latex')
      title('$$p_A(a_1,a_2)$$','interpreter','latex')
      
      figure; plot(XevalB{1},pB);
      xlabel('$$b_1$$','interpreter','latex')
      ylabel('$$p_B(b_1)$$','interpreter','latex')
    else
      ylabel('$$p_B(b_2)$$','interpreter','latex')
      zlabel('$$p_B(b_3)$$','interpreter','latex')
      
      figure; plot(XevalA{1},pA);
      xlabel('$$a_1$$','interpreter','latex')
      ylabel('$$p_A(a_1)$$','interpreter','latex')
      
      figure; imagesc(XevalB{1},XevalB{2},pB);
      xlabel('$$b_1$$','interpreter','latex')
      ylabel('$$b_2$$','interpreter','latex')
      title('$$p_B(b_1,b_2)$$','interpreter','latex')
    end
  else
    
  end
end


