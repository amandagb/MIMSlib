function [pdfkde,p_i,clk,animals,pdfhist,Imasks,expinfo] = EB_brainimg_buildpdfs(varargin);
%% function: EB_brainimg_buildpdfs
% [pdfkde,p_i,clk] = EB_brainimg_buildpdfs(animals,pdfhist,Imasks,expinfo);
% Description: Reads in data generated from EB_brainimg_processing function
% and then constructs class pdfs using kde as well as rgb intensity pdfs
% INPUTS ----------------------------------------------------------------
% varargin{1}: date string in the form of 'YYYYMMDD' [DEFAULT = '20141030']
% varargin{1:3}: outputs from EB_brainimg_processing function
%   {1}: animals, {2}: pdfhist, {3}: Imasks, {4}: exp info
% varargin - 'PropertyName','PropertyValue'
%   *• 'nclasses': scalar indicating whether to use 2 classes (pathology,
%   normal tissue), 3 (pathology, normal, background), 4 (blood, EB,
%   complex, normal), or 5 (all in 4 + background) [DEFAULT = 4]
%   * Not yet included in code
%   • 'kdebw': scalar indicating the desired bandwidth of the kernel. In
%   the case of symmetric normal in all dimensions, this number represents
%   the VARIANCE in a single dimension
% 
% OUTPUTS ---------------------------------------------------------------
% pdfkde:     structure with fields containing 256x256x256 doubles which
% are each pdfs for that class => each struct = p(i|c)
% p_i:        structure with fields containing 256x256x256 doubles which
% are each pdfs for that class => each struct = p(i|fgnd only) or p(i)
% 
%  Date           Author              E-mail                      Version
%   8 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%     Isolated as a separate function from EB_brainimg_classifyimgs

close all;
loaddir;

if length(varargin) == 1
  startvin = 2;
  data_date = varargin{1};
  EBdpath = strcat(dbEBexppath,data_date,'_EBclasshistmed');
  load(EBdpath,'animals','pdfhist','Imasks','expinfo');
elseif length(varargin) == 4
  startvin = 5;
  animals = varargin{1};
  pdfhist = varargin{2};
  Imasks = varargin{3};
  expinfo = varargin{4};
else
  startvin = 1;
  data_date = '20141030';%'20141026';%'20141016';
  EBdpath = strcat(dbEBexppath,data_date,'_EBclasshistmed');
  load(EBdpath,'animals','pdfhist','Imasks','expinfo');
end
data_date = expinfo.data_date;
EBprocessingversion = expinfo.version;

% if exist(EBdpath)
  
% else 
%   disp(sprintf('There is no file associated with the date string indiciated (%s)',data_date));
%   return;
% end
loaddir;
EBfldr = strcat(fpath,'EB_brainimg');

PropertyNames = varargin(startvin:2:length(varargin));
PropertyVal = varargin(startvin:2:length(varargin));

if strmatch('kdebw',PropertyNames)
  norm_var = PropertyVal{strmatch('kdebw',PropertyNames)};
else
  norm_var = 9;
end

%=========================================================================
% Initialize folder & File Variables
%=========================================================================
extraclasses = pdfhist.otherclasses; %{'bkgnd','rgb','frgb'};
class_str = pdfhist.classes; %{'bld','EB','hemo','tsu'};
full_class_str = pdfhist.fullclasses; %{'Blood','EB','Contusion','Tissue'};
class_color = {'r','b','k','g'};

animalnumvec = animals.id;%[animals.id(1)];%
plotimgs = 1;

%=========================================================================
% Split data into relevant subgroups: KDE is not conducted for prob. rgb
% intensity. This will llikely need to be changed once only training data
% is used to construct these densities unless they are assumed to be a
% "state of nature" type distribution 
%=========================================================================
% Create smaller data: background class, all rgb class, rgb in fgnd only
bkgndhist = getfield(pdfhist,'bkgnd');
rgbhist = double(getfield(pdfhist,'rgb'));
frgbhist = double(getfield(pdfhist,'frgb'));
p_i.bkgnd = bkgndhist;  % OUTPUT
p_i.all = rgbhist(:)./sum(rgbhist(:));      % OUTPUT
p_i.fonly = frgbhist(:)./sum(frgbhist(:));  % OUTPUT

use4class = pdfhist;
use4class = rmfield(use4class,extraclasses);
use4class = rmfield(use4class,class_str);
classhist = rmfield(pdfhist,extraclasses);%,'frgb'});
classhist = rmfield(classhist,fieldnames(use4class));%,'frgb'});
fnames = fieldnames(classhist);
nclasses = length(fnames);
show_perc_data = 100; % percentage of data that should be displayed on the plot
pdfkde = cell(1,nclasses);
w = 1; 
for c = 1:nclasses;
  clk(w,:) = clock; % w = 1
  w = w + 1;
  data = getfield(pdfhist,fnames{c});
  nhistbins = length(data(:));
  npnts = double(sum(data(:)));
  nnon0bins = double(sum(data(:) ~= 0));
  disp(sprintf('%s: %d data points, %d nonzero bins of %d (%1.2f %%)',fnames{c},npnts,nnon0bins,nhistbins,nnon0bins/nhistbins*100))
  % One VERY IMPORTANT point is that U(1) will almost definitely be 0. This
  % should be disregarded for plotting and other counting purposes.
  % Therefore, ic == 1 should also be disregarded, but does tell how many
  % unfilled bins (or completely absent colors) there are in the rgb histogram
  [U,ia,ic] = unique(data);
  nU = length(U);
  U = double(U);
  Ufrac = U./npnts;
  Ucs = cumsum(U)./sum(U);
  COth = 1 - show_perc_data/100; %cut off threshold -- this is fraction of data that will be cut off
  Uismall = find(Ucs < COth);
  if isempty(Uismall)
    uicutoffth = 0;
  else
    uicutoffth = Ufrac(Uismall(end)); % cutoff threshold --- this plots only bins that have great than cutoff thresh fraction of the total data points in it
  end
  uicutoff = find(U > npnts*uicutoffth);
  clk(w,:) = clock; % w = 2
  t_defbins(c) = etime(clk(w,:),clk(w-1,:));
  w = w + 1;
  
  %=========================================================================
  % Visit each of the non-zero bins and drop a kernel in the bin. weigh by
  % the count of the data in the bin
  %=========================================================================
  %{.
  % ------------------------------------------------------------------------
  % • Parameterize the KDE kernel K -- most parameters identified from Hansen
  % 2009 "Lecture Notes on Nonparametrics"
  %   nu: kernel order
  %   RK: kernel roughness
  %   kap2: second moment of the kernel
  %     - Parameters are given on table 1, page 4 of the cited document
  % ------------------------------------------------------------------------
  d = 3;
  Kdim{1} = 1:d;
  Kstr{1} = 'norm';
  nu = 2; % kernel order
  Kcell = cell(d,1);
  bw = zeros(1,d);
  Ksupp = nan(2,d); % support of the kernel
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
    % --> Rather than computing an automatically optimal bw using the
    % data's standard deviation, the standard deviation is fixed to be 15
    %     if isempty(Kbw{kd})
    %       sig_hat = std(rgbpnts);
    %       CK = 2*((pi^(1/2)*factorial(nu)^3*RK)/(2*nu*factorial(2*nu)*kap2^2))^(1/(2*nu+1));
    %       bw = sig_hat.*CK.*npnts^(-1/(2*nu+1));
    %     else
    %       bwkd = Kbw{kd};
    %       if length(bwkd) < length(Kdim{kd})
    %         bw(Kdim{kd}) = repmat(bwkd,1,length(Kdim{kd}));
    %       else
    %         bw(Kdim{kd}) = bwkd;
    %       end
    %     end
  end
  
  Ksupp = Ksupp.*repmat(bw,2,1);
  %   clk(w,:) = clock; % w = 5
  %   Kformtime = etime(clk(w,:),clk(w-1,:));
  %   w = w + 1;
  %}
  
  %---- Define kernel BW
  sig_hat = repmat(15,d,1); %sig_hat = std(rgbpnts);
  CK = 2*((pi^(1/2)*factorial(nu)^3*RK)/(2*nu*factorial(2*nu)*kap2^2))^(1/(2*nu+1));
  bw = sig_hat.*CK.*npnts^(-1/(2*nu+1));
  
  %--- Automatically determine kernel support required by cutting off all
  %values that would result in something < KCOlim AFTER being scaled
  % SEE LATEX notes 20141101mtng notes for details. Essentially, the
  % problem is solved easily since we assume no covariance and equal
  % variance. - a symmetric Gaussian kernel
  maxbincnt = max(U);
  KCOlim = 1e-6;
  bw_var = repmat(norm_var,1,d - size(norm_var,2) + 1);
  S = diag(bw_var);
  Sinv = S^(-1);
  norm_dconst = sqrt((2*pi)^d * det(S));
  LHSconst = -2*log((KCOlim*norm_dconst)/maxbincnt);
  % Here, since the gaussian is symmetric in all directions, we can assume
  % that varlim will be equal for all axes
  varlim = sqrt(LHSconst/Sinv(1,1)); 
  varlimint = ceil(varlim);
  
  rgbrng = [-varlimint:varlimint];
  [R,G,B] = meshgrid(rgbrng,rgbrng,rgbrng);
  Ktemp = mvnpdf([R(:),G(:),B(:)],zeros(1,d),S);
  Ktemp = reshape(Ktemp,[size(R,1),size(R,2),size(R,3)]);
  
  clk(w,:) = clock; % w = 3
  w = w + 1;
  cpdf = zeros(256+2*varlimint,256+2*varlimint,256+2*varlimint);
  disp(sprintf('Eval %d bin cnts, Max bin cnts = %d, kde range = +/- %d',...
    length(uicutoff), maxbincnt, varlimint))
  for i = uicutoff(:)'
    binscale = U(i);  
    icxUi = find(ic == i); % identifies the indices which correspond to the bin scale
    [ri,gi,bi] = ind2sub([256,256,256],icxUi); % indices mapped into RGB space
    evali = [ri,gi,bi]; % creates a # indices x 3 vector
    newi = evali + varlimint; % offset associated with indices since the 256x256x256 matrix was padded on all sides
    for p = 1:size(newi,1) %
      pntrng = repmat(newi(p,:),length(rgbrng),1) + repmat(rgbrng',1,d);
      cpdf(pntrng(:,1),pntrng(:,2),pntrng(:,3)) =...
        cpdf(pntrng(:,1),pntrng(:,2),pntrng(:,3)) + Ktemp.*binscale;
    end
  end
  
  cpdf = cpdf(varlimint+(1:256),varlimint+(1:256),varlimint+(1:256));
  cpdf = cpdf./sum(cpdf(:));
  if ~exist('pdfkde')
    pdfkde = struct(fnames{c},cpdf);
  else
    pdfkde = setfield(pdfkde,fnames{c},cpdf);
  end

  clk(w,:) = clock; % w = 4
  tclass(c) = etime(clk(w,:),clk(w-1,:));
  disp(tclass(c))
  w = w + 1;
end

end