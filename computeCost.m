function [J,totpxls] = computeCost(F,M,coststr,varargin)
%% Script: computeCost
% Description:
% Example:
% Required Functions: twoimgMIkde
% INPUTS ----------------------------------------------------------------
% F:  M x N x df fixed image (typically taken to be double precision)
% M:  M x N x dm moving image
% coststr: string indicating the cost to compute
%   • 'mi'
%   • 'mine'
%   • 'lp'
%   • 'cc'
%   • 'div'
%   • 'tv'
%   • 'max1'
%   • 'HL'
% varargin - 'PropertyName','PropertyValue'
%   • 'mu0':      1x6 vector with initial transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   • 'pval':     scalar indicating the p-value of the ell-p norms [DEFAULT = 1]
%   • 'truncate': scalar indicating maximum value of the cost function
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%   8 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1
%     Based on v4 of simanneal_MI
%  13 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1.1
%     Added scaling of cost function

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

[Mm,Mn,Md] = size(M);
[Fm,Fn,Fd] = size(F);
%.....................................................................
%   varargin
%.....................................................................

if strmatch('Fmask',PropertyNames)
  Fmask = PropertyVal{strmatch('Fmask',PropertyNames)};
else
  Fmask = true(Fm,Fn);
end

if strmatch('Mmask',PropertyNames)
  Mmask = PropertyVal{strmatch('Mmask',PropertyNames)};
else
  Mmask = true(Mm,Mn);
end

if strmatch('pval',PropertyNames)
  pval = PropertyVal{strmatch('pval',PropertyNames)};
else
  pval = 1; %tx, ty, theta, sx, sy, sk
end

if strmatch('mu',PropertyNames)
  mu = PropertyVal{strmatch('mu',PropertyNames)};
else
  mu = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk
end

if strmatch('MImethod',PropertyNames)
  MImethod = PropertyVal{strmatch('MImethod',PropertyNames)};
else
  MImethod = 'kde';
end

if strmatch('sumdim',PropertyNames)
  sumdim = PropertyVal{strmatch('sumdim',PropertyNames)};
else
  sumdim = 0;
end
%.....................................................................
%   data preparation
%.....................................................................
[Mm,Mn,Md] = size(M);
[Fm,Fn,Fd] = size(F);
Fvec = reshape(F(:),[Fm*Fn,Fd]);

Fnan = F;
Fnan(~repmat(Fmask,[1,1,Fd])) = nan;

Mnan = M;
Mnan(~repmat(Mmask,[1,1,Md])) = nan;

TxM = transform_image(M,mu,'extrapval',M(1));
TxMask = transform_image(double(Mmask),mu,'extrapval',double(Mmask(1)));
usei = and(Fmask,TxMask);
totpxls = sum(usei(:)) - sum(sum(isnan(sum(cat(3,Fnan,Mnan),3))));
TxMvec = reshape(TxM(:),[Mm*Mn,Md]);
if Fd > Md
  dFM = Fvec - [TxMvec,zeros(Fm*Fn,Fd - Md)];
elseif Md > Fd
  dFM = [Fvec,zeros(Fm*Fn,Md - Fd)] - TxMvec;
else
  dFM = Fvec - TxMvec;
end

switch coststr
  case 'mi'
    switch MImethod
      case 'mattes'
        
      case 'kde'
        if sumdim
          if Fd > 1
            for i = 1:Fd
              if nansum(nansum(Fnan(:,:,i)))
                MI(i) = twoimgMIkde(Fnan(:,:,i),TxM,varargin{:});%'nbins',nbins);
              else
                MI(i) = 0;
              end
            end
          else
            [MI,outinfo] = twoimgMIkde(Fnan,TxM,varargin{:});%'nbins',nbins);
          end
          J = -sum(MI);
        else
          [MI,outinfo] = twoimgMIkde(Fnan(:,:,1:Fd),TxM,varargin{:});%,'mu',mu);%'nbins',nbins);
          totpxls = outinfo.nonNANpxls;
          J = -MI;
        end
    end
  case 'lp'
    %dvec = Fvec(usei(:),:) - TxMvec(usei(:),:);
    %J = norm(dvec(:),pval); --- use of matlab native function could be more efficient
    dFM = sum(abs(dFM).^pval,2);
    Jimg = reshape(dFM,[Fm,Fn]);
    J = sum(Jimg(usei));
  case 'mine' % Albanese, minepy-1.0.0, Bioinformatics (2013) see: minepy.sourceforge.net (based on Reshef 2011 MIC metric paper)
    mineout = mine(Fvec(usei(:),:)',TxMvec(usei(:),:)');
    J = mineout.mic;
  case 'lcca' % Heinrich, MICCAI 2014, see http://www.mpheinrich.de/software.html
    Jimg = LCCA_2d(F, TxM, 5, 1e-3);
    J = sum(Jimg(usei));
  case 'cca'
    [~,~,Jimg,~,~,~] = canoncorr(Fvec(usei(:),:),TxMvec(usei(:),:));
    J = sum(Jimg);
  case 'div' % J(x,x0) = sum_i x_i log(x_i/x_0i) + sum_i (x_i - x_0i), normalize: sum_i x_i = 1, uniform prior: x_0i = 1/N => J(x) = sum_i x_i log x_i
    %x = abs(dFM(usei))./sum(abs(dFM(usei)));
    J = entropy(dFM(usei))*log(2); % in NATs
  case 'tv'
    dFM = sum(abs(dFM).^pval,2);
    Jimg = reshape(dFM,[Fm,Fn]);
    J = nanmax(dFM) - nanmin(dFM);
  case 'max1'
    dFM = sum(abs(dFM).^pval,2);
    Jimg = reshape(dFM./(1+dFM),[Fm,Fn]);
    J = nansum(Jimg(usei));
  case 'HL'
    dFM = sum(log(abs(dFM).^pval + 1),2);
    Jimg = reshape(dFM,[Fm,Fn]);
    J = sum(Jimg(usei));
end
end