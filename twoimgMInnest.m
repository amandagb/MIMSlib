function [phi_MI,varargout] = twoimgMInnest(A,B,varargin)
%% Script: [phi_MI,varargout] = twoimgMInnest(A,B,varargin)
% Description: Nearest-neighbor pdf/functional entropy estimation
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
%   • 'mc':    # Monte Carlo runs
%   • 'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations [DEFAULT = 50]
%   • 'logbase': string indicating whether the log base should be '2' or 'e'
%   [DEFUALT = 'e']
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
%   Multiple kernels can be chosen by letting the varaible have multiple
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
%  13 Nov  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%     Base script derrived from twoimgMIkde (v0.2)

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

if strmatch('mc',PropertyNames)
  mcruns = PropertyVal{strmatch('mc',PropertyNames)};
else
  mcruns = 100;
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

c = zeros(mcruns + 20,6);
% ------------------------------------------------------------------------
% • Turns image A from an [m x n x d] image into an [mn x d] vector
% • Defines the bin centers
% ------------------------------------------------------------------------
w = 1; c(w,:) = clock; % w = 1
w = w + 1;

if isempty(Abinlim)
  minA = min(Aflat);
  maxA = max(Aflat);
  Abinlim(:,1) = minA';
  Abinlim(:,2) = maxA';
else
  minA = Abinlim(:,1)';
  maxA = Abinlim(:,2)';
end

for k = 1:dA
  Aflat(Aflat(:,k) < minA(k),k) = nan;
  Aflat(Aflat(:,k) > maxA(k),k) = nan;
end

c(w,:) = clock; % w = 3
Aformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Turns image B from an [m x n x d] image into an [mn x d] vector
% • Transforms B according to indicate transformation parameters mu
% • Defines the bin centers
% ------------------------------------------------------------------------
TxB = transform_image(B,mu,varargin{:});%,'interpmethod','nearest'); % added property 'outputimsz' 2013 April 9
TxBflat = reshape(TxB(:),MrB*NcB,dB);
if isempty(Bbinlim)
  maxB = max(Bflat);
  minB = min(Bflat);
  Bbinlim(:,1) = minB';
  Bbinlim(:,2) = maxB';
else
  minB = Bbinlim(:,1)';
  maxB = Bbinlim(:,2)';
end

for k = 1:dB
  TxBflat(TxBflat(:,k) > maxB(k),k) = nan;%maxB(k);
  TxBflat(TxBflat(:,k) < minB(k),k) = nan;%minB(k);
end

c(w,:) = clock; % w = 3
Bformtime = etime(c(w,:),c(w-1,:));
w = w + 1;

% ------------------------------------------------------------------------
% • Defines joint sample points Xi
% • Defines function evaluation points Xeval
% ------------------------------------------------------------------------
Z = [Aflat, TxBflat];
useXi = find(sum(isnan(Z),2) == 0);
Z = Z(useXi,:);
Aflat = Aflat(useXi,:);
TxBflat = TxBflat(useXi,:);
npnts = length(useXi);
m_pdfest = round(npnts/2);
n_fctnest = npnts - m_pdfest;

% ------------------------------------------------------------------------
% • Parameterize the KDE kernel K
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
      %Ksupp(:,Kdim{kd}) = [-1;1];
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
    sig_hat = std(Z);
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
Hest = zeros(mcruns,3);
Htimes = zeros(mcruns,3);
for mc = 1:mcruns
  sampind = randsample(npnts,npnts); % Should these be generated randomly for each random variable or should the same indices be used for estimating pdf & fctnl for all random variables?
  pdfi = sampind(1:m_pdfest);
  fctnli = sampind(m_pdfest+1:end);
  for rv = 1:3
    switch rv
      case 1
        RV = Aflat;
        jointd = 1:dA;
        %bw = std(Aflat);
        %muRV = mu0(1);
        %sigRV = sig2x;
      case 2
        RV = TxBflat;
        jointd = (dA+1):(dA+dB);
        %bw = std(TxBflat);
        %muRV = mu0(1);
        %sigRV = sig2y;
      case 3
        RV = Z;
        jointd = 1:(dA+dB);
        %bw = std([X,Y]);
        %muRV = mu0(1);
        %sigRV = Sig1;
    end
    tic;
    X_pdfest = RV(pdfi,:);
    X_fctnest = RV(fctnli,:);
    
    outsidearg = 0;
    for i = 1:n_fctnest
      %insidearg = 0;
      Xi = X_fctnest(i,:);
      dXi = abs(repmat(Xi,m_pdfest,1) - X_pdfest); % order doesn't matter since all kernels are symmetric
      
      for kd = 1:size(Kinfo,1)
        datad = Kdim{kd};
        memi = ismember(jointd,datad);
        filli = find(memi == 1);
        Kout(:,filli) = (Kcell{datad(1)}(dXi(:,filli)...
          ./repmat(bw(filli),m_pdfest,1)))./(repmat(bw(filli),m_pdfest,1));
      end
      
      outsidearg = outsidearg + (-log(sum(prod(Kout,2))/m_pdfest))/n_fctnest;
      %outsideC = outsideC + (-log(mvnpdf(Xi,muRV,sigRV)))/n_fctnest;
    end
    Hest(mc,rv) = outsidearg;
    %C(mc,samp,rv) = outsideC;
    c(w,:) = clock; % w = 6
    Htimes(mc,rv) = etime(c(w,:),c(w-1,:));
    w = w + 1;
  end
end

% ------------------------------------------------------------------------
% • Actual computation of MI by summing the MI arguments generated before
% ------------------------------------------------------------------------
phi_MI = Hest(:,1) +  Hest(:,2) - Hest(:,3);

c(w,:) = clock; % w = 13
MIargsumtime = etime(c(w,:),c(w-1,:));

tot_time = etime(c(w,:),c(1,:));

all_times.Aformtime = Aformtime;
all_times.Bformtime = Bformtime;
% all_times.Xevaltime = Xevaltime;
all_times.Kformtime = Kformtime;
all_times.HAtime = Htimes(:,1);
all_times.HBtime = Htimes(:,2);
all_times.HABtime = Htimes(:,3);
all_times.MIargsumtime = MIargsumtime;
all_times.tot_time = tot_time;

dist.Abinlim = Abinlim;
dist.Bbinlim = Bbinlim;
dist.HA = Hest(:,1);
dist.HB = Hest(:,2);
dist.HAB = Hest(:,3);
dist.mcruns = mcruns;
dist.totpnts = npnts;
dist.pdfpnts = m_pdfest;
dist.fctnpnts = n_fctnest;
dist.time = all_times;
dist.kerninfo.bw = bw;
dist.kerninfo.str = Kstr;

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


 