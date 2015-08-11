function [phi_MI,varargout] = twoimgMImd(A,B,varargin)
%% Script: [phi_MI] = twoimgMImd(A,B,varargin)
% Description: The following implementation is derived from
%   [Mattes et al.(2003)] D. Mattes, D. Haynor, H. Vesselle, T. Lewellen,
%   and W. Eubank. PET-CT image registration in the chest using free-form
%   deformations. Medical Imaging, IEEE Transactions on, 22(1):120 –128,
%   Jan 2003. ISSN 0278-0062. doi: 10.1109/ TMI.2003.809072.
%   Multi-dimensional version of twoimgMImd using multidimensional
%   b-splines
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% A:  MrA x NcA Image with dA channels (fixed/reference image). May also be
% a MrA·NcA x dA vector. if NcA is < 10, it will automatically be treated
% as the latter.
% B:  MrB x NcB Image with dB channels (moving/test image) -- see A
% varargin - 'PropertyName','PropertyValue'
%   'mu':    1x6 vector with transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations [DEFAULT = 50]
%   'logbase': string indicating whether the log base should be '2' or 'e'
%   [DEFUALT = '2']
%   'pA': 1 x nbins vector which approximates the pdf of the fixed image, A
%   'plot': binary indicating whether to plot (1) or not (0) [DEFAULT = 0]
%   'axh': axis handle which should be used for plotting
%   'Abinlim': dA x 2 vector indicating the [min,max] bin edge values for
%   each channel
%   'Bbinlim': dB x 2 vector indicating the [min,max] bin edge values for
%   each channel
% OUTPUTS ---------------------------------------------------------------
% phi_MI: Mutual information between fixed image A and moving image B
% varargout
%   {1} structure with fields which contain the distributions for A (pA), B
%   (pB) and the joint (pAB)
%       pA: distribution of A
%       Alims: bin centers where values of A range from min(A) to max(A)
%       and bin width can be computed as (Alims(1)-min(A))*2
%       pB: distribution of B
%       Blims: bin centers where values of A range from min(A) to max(A)
%       and bin width can be computed as (Alims(1)-min(A))*2
%       pAB: joint distribution
%
%
%  Date           Author            E-mail                      Version
%  15 Aug 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%     Initial draft/debugging script (not a function)
%  26 Aug 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     Exhaustively computes pAB for each [k,l] set (for every point in pAB)
%  26 Aug 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Computes pAB only for existing ABjoint values
%  29 Aug 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Added clocking and various methods for improving computation time
%     > Vectorization of pdf's was attempted in v0.1, see OldScripts
%     > v0.2: Attempted to exclude sparsely populated bins from the
%     calculation
%     > [18 Sept] v3.3: added additional variable input which allows user
%     to specify the range of the bins
%     > [20 Sept] v3.4: Fixed the binning

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
  logbase = '2';
end

if strmatch('nbins',PropertyNames)
  nbinsA = PropertyVal{strmatch('nbins',PropertyNames)}.*ones(1,dA);
  nbinsB = PropertyVal{strmatch('nbins',PropertyNames)}.*ones(1,dB);
else
  nbinsA = 50.*ones(1,dA);
  nbinsB = 50.*ones(1,dB);
end

if strmatch('nbinsA',PropertyNames)
  nbinsA = PropertyVal{strmatch('nbinsA',PropertyNames)};
end

if length(nbinsA) < dA
  Abins = nbinsA.*ones(1,dA);
else
  Abins = nbinsA;
end

if strmatch('nbinsB',PropertyNames)
  nbinsB = PropertyVal{strmatch('nbinsB',PropertyNames)};
end

if length(nbinsB) < dB
  Bbins = nbinsB.*ones(1,dB);
else
  Bbins = nbinsB;
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

% ------------------------------------------------------------------------
% • Turns image A from an [m x n x d] image into an [mn x d] vector
% • Implements the zero order B-spline definition in the for loop to
% represent the vector A image as bin indices rather than original
% intensity values
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
bin_rngA = Arng./Abins;
Aind_norm = (Aflat - repmat(minA,size(Aflat,1),1))./repmat(Arng,size(Aflat,1),1); %repmat(bin_rngA,size(Aflat,1),1);
%%% ANOTHER OPTION HERE IS TO SIMPLY EXCLUDE THESE INDICES WHICH FALL
%%% OUTSIDE THE RANGE FROM THE FINAL COMPUTATION
Aind_norm(Aind_norm < 0) = nan;
Aind_norm(Aind_norm > 1) = nan;
Aind_float = Aind_norm.*repmat((Abins -1),size(Aflat,1),1);
Aind_int = ceil(Aind_float);
Aind_int(Aind_int == 0) = 1;
Alims = cell(1,dA);
for k = 1:dA
  Alims{k} = [0.5:1:(Abins(k)-0.5)].*bin_rngA(k) + minA(k);
end

c(w,:) = clock; % w = 2
Aformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Turns image B from an [m x n x d] image into an [mn x d] vector
% • Implements the zero order B-spline definition in the for loop to
% represent the vector B image as bin indices rather than original
% intensity values -- these will be used to guide the relevant bins for the
% 3rd order B-spline
% ------------------------------------------------------------------------
TxB = transform_image(B,mu,'extrapval',nan);%,'outputimsz',size(A));%,'interpmethod','nearest'); % added property 'outputimsz' 2013 April 9
TxBflat = reshape(TxB(:),MrB*NcB,dB);
if isempty(Bbinlim)
  maxB = max(Bflat);
  minB = min(Bflat);
else
  minB = Bbinlim(:,1)';
  maxB = Bbinlim(:,2)';
end
Brng = maxB - minB;
bin_rngB = Brng./Bbins;
Blims = cell(1,dB);
for k = 1:dB
  TxBflat(TxBflat(:,k) > maxB(k),k) = nan;%maxB(k);
  TxBflat(TxBflat(:,k) < minB(k),k) = nan;%minB(k);
  Blims{k} = [0.5:1:(Bbins(k)-0.5)].*bin_rngB(k) + minB(k);
end
if isempty(Bbinlim)
  maxB = max(TxBflat);
  minB = min(TxBflat);
end
Bind_float = (TxBflat - repmat(minB,size(TxBflat,1),1))./repmat(bin_rngB,size(TxBflat,1),1);
Bind_int = ceil(Bind_float);
Bind_int(Bind_int == 0) = 1;

c(w,:) = clock; % w = 3
Bformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Identify unique AB-bin pairs
% ------------------------------------------------------------------------
binsvec = [Abins,Bbins];
pAB = zeros(binsvec); %zeros(prod(binsvec),1);
AintBfloat = [Aind_int,Bind_float];
AintBfloat = AintBfloat(isnan(sum(AintBfloat,2))== 0,:);
npnts = size(AintBfloat,1);

Ajoints = unique(Aind_int,'rows');
ABjoints = unique([Aind_int,Bind_int],'rows');
ABjoints = ABjoints(isnan(sum(ABjoints,2))== 0,:);

c(w,:) = clock; % w = 4
unqABpairtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Defines unique AB-bin pairs that need to be evaluated by defining
% relevant neighborhood around B's indices which can potentially be filled
% by kernel
% ------------------------------------------------------------------------
repind = repmat([1:size(ABjoints,1)],5^dB,1);
ABj = ABjoints(repind(:),:);
j = 1;
Bpossible = combvec([-2:2]);
while j < (dB)
  Bpossible = combvec(Bpossible,[-2:2]);
  j = j+1;
end
Bpossible = Bpossible';
ABj = ABj + [zeros(size(ABj,1),dA), repmat(Bpossible,size(ABjoints,1),1)];
ABjoints = unique(ABj,'rows');
ABjoints = ABjoints(sum(sign(ABjoints-1)<0,2) == 0,:);
ABjoints = ABjoints(sum(ABjoints < repmat(binsvec+1,size(ABjoints,1),1),2) == (dA+dB),:);

c(w,:) = clock; % w = 5
totABpairtime = etime(c(w,:),c(w-2,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute pAB by considering only existing ABjoint values
% ------------------------------------------------------------------------
ABjlinind = ones(size(ABjoints,1),1);

for j = 1:(dA+dB)
  ABjlinind = (ABjoints(:,j)-1).*binsvec(j)^(j-1) + ABjlinind;
end

for kv = 1:size(Ajoints,1)
  k = Ajoints(kv,:);
  krows = find(sum(repmat(k,npnts,1) == AintBfloat(:,1:dA),2)==dA);
  relpnts = length(krows);
  Brel = AintBfloat(krows,dA+1:end);
  ABji = find(ismember(ABjoints(:,1:dA),k,'rows'));
  
  %   skip = zeros(length(ABji),1);%   zerol = zeros(1,dB);
  %   n = 0;
  for i = ABji'
    %     n = n + 1;
    l = ABjoints(i,dA+1:end);
    %     if ~skip(n)
    Barg = abs(repmat(l-1,relpnts,1) - Brel);
    Bargvalid = Barg > 2;
    non0rows = find(sum(Bargvalid,2) == 0);
    if ~isempty(non0rows)
      Barg = Barg(non0rows,:);
      Afunvec = ones(size(Barg,1),1);
      
      Bfunvec = zeros(size(Barg));
      Bfunvec(Barg <= 1) = 1/6.*(4 - 6.*Barg(Barg <= 1).^2 + 3.*Barg(Barg <= 1).^3);
      Bfunvec(Barg < 2 & Barg > 1) = 1/6.*(2-Barg(Barg < 2 & Barg > 1)).^3;
      pAB(ABjlinind(i)) = sum(prod(Bfunvec,2).*Afunvec);
      %       elseif any(sum(Bargvalid,1) == relpnts)
      %         badlcol = find(sum(Bargvalid,1) == relpnts);
      %         skip = ABjoints(ABji,badlcol(1) + dA) == l(badlcol(1));
    end
    %     end
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
Asumdim = (dA+dB):-1:(dA+1);
pA = pAB;
for i = Asumdim
  pA = sum(pA,i);
end

c(w,:) = clock; % w = 8
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

c(w,:) = clock; % w = 9
pBcalctime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Compute denominator of MI log function
% ------------------------------------------------------------------------
pkpl_product = pA(:)*pB(:)';%reshape(pA(:)*pB(:)',binsvec);
% Another opportunity to speed up computation is to look only at the pA, pB
% product pairs which will be used in MIarg. That is, only compute
% length(non0pAB) products which correspond to the nonzero indices in pAB

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

dist.pABnon0 = pAB(non0pAB);
dist.logarg = (pAB(non0pAB))./(pkpl_product(non0pAB));
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
    ax1 = subplot('position',[0.2,0.2,0.7,0.7]); imagesc(pAB);
    set(gca,'xticklabel','','yticklabel','')
    subplot('position',[0.1,0.2,0.1,0.7]); plot(pA,1:Abins,'marker','.'); axis ij;
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
    subplot('position',[0.2,0.1,0.7,0.1]); plot(1:Bbins,pB,'marker','.');
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
  elseif dA + dB == 3
    hold on;
    cc = custom_colormaps('jet');
    colind = floor(pAB./max(pAB(:)).*size(cc,1)) + 1;
    colind(find(colind == size(cc,1) + 1)) = size(cc,1);
    ppbthresh = 1; % points per bin threshold
    i = 1;
    pall = combvec(1:binsvec(i));
    while i < (dA+dB)
      pall = combvec(pall,1:binsvec(i+1));
      i = i+1;
    end
    
    for i = 1:length(pall)
      
      linind = 1;
      for j = 1:(dA+dB)
        linind = (pall(j,i)-1)*binsvec(j)^(j-1) + linind;
      end
      if pAB(pall(1,i),pall(2,i),pall(3,i)) > ppbthresh/numel(pAB)
        plot3(pall(1,i),pall(2,i),pall(3,i),'s','markersize',10,...
          'markerfacecolor',cc(colind(linind),:),...
          'markeredgecolor',cc(colind(linind),:));
      end
    end
    xlabel('$$p_A(x_1)$$','interpreter','latex')
    %     Aticks = str2num(get(gca,'xtickLabel'));
    %     ntick = length(Aticks);
    %     ntick = 11;
    %     tickdist = round((max(Aflat) - min(Aflat))./(ntick-1));
    %     Aticklab = floor(min(Aflat)):tickdist:ceil(max(Aflat));
    %     Aticklab(1) = ceil(min(Aflat(:,1))*10)/10;
    %     Aticklab(end) = floor(max(Aflat)*10)/10;
    %     tickpos = (Aticklab-min(Aflat))./bin_rngA; % index image which reference which bin the pixel belongs to
    %     set(gca,'xtick',tickpos,'xticklabel',Aticklab);
    
    if dA == 2
      ylabel('$$p_A(x_2)$$','interpreter','latex')
      zlabel('$$p_B(x_3)$$','interpreter','latex')
      
      figure; imagesc(pA);
      xlabel('$$x_1$$','interpreter','latex')
      ylabel('$$x_2$$','interpreter','latex')
      title('$$p_A(x_1,x_2)$$','interpreter','latex')
      
      figure; plot(pB);
      xlabel('$$x_1$$','interpreter','latex')
      ylabel('$$p_B(x_1)$$','interpreter','latex')
    else
      ylabel('$$p_B(x_2)$$','interpreter','latex')
      zlabel('$$p_B(x_3)$$','interpreter','latex')
      
      figure; plot(pA);
      xlabel('$$x_1$$','interpreter','latex')
      ylabel('$$p_A(x_1)$$','interpreter','latex')
      
      figure; imagesc(pB);
      xlabel('$$x_1$$','interpreter','latex')
      ylabel('$$x_2$$','interpreter','latex')
      title('$$p_B(x_1,x_2)$$','interpreter','latex')
      %       Bticks = str2num(get(gca,'xtickLabel'));
      %       ntick = length(Bticks);
      %       ntick = 11;
      %       tickdist = round((max(TxB(:)) - min(TxB(:)))/(ntick-1));
      %       Bticklab = floor(min(TxB(:))):tickdist:ceil(max(TxB(:)));
      %       Bticklab(1) = ceil(min(TxB(:))*10)/10;
      %       Bticklab(end) = floor(max(TxB(:))*10)/10;
      %       tickpos = (Bticklab-min(TxB(:)))./bin_rngB; % index image which reference which bin the pixel belongs to
      %       set(gca,'xtick',tickpos,'xticklabel',Bticklab);
      %       %useBticks = Bticks(Bticks > 0);
      %       %set(gca,'xtick',useBticks,'xticklabel',Blims(useBticks));
      %       xlabel('Moving image intensities, $$y$$','interpreter','latex')
      %       ylabel('$$p(y)$$','interpreter','latex')
      
    end
  else
    
  end
end
%end