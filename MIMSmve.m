function [dout,mveparams,varargout] = MIMSmve(data,useElstr,varargin);
%% Script:  data_struct = MIMSmve;
% Description:
%
% Example:
%   >> [dprocMF,mveparMF] = MIMSmve('DATASTR',useElstr,'labelstr','HIPP',...
%         'medfilter',0,'repback',1);
%
% Required Function:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% elem_ana_range: 1 x <# elements> matrix which corresponds to the field
%                 number in the data structure. User may also specify
%   This variable can also be a cell of strings where the string specifies
%   the isotope the user would like to analyze
% varargin - 'PropertyName','PropertyValue'
%   • 'hdrtxt': cell containing the appropriate data hdrtxt variable
%   • 'labelstr': string indicating the corresponding label map to isolate
%   and plot
%   • 'repback': logical indicating whether to replace the background with
%   a single scalar value according to tissue mask (1) or not (0) [DEFAULT = 1]
%   • 'rotdat': logical indicating whether to rotate data (1) or not (0) [DEFAULT = 0]
%   • 'plotimgs': logical indicating whether to plot images (1) or not (0)
%   [DEFUALT = 0]
%   • 'saveimgs': logical indicating whether to save images (1) or not (0)
%   [DEFUALT = 0]
%   • 'outtag': output data tag to indicate the type of transformed data
%   that should be output
%       > 'raw': Do no processing (unless 'repback' is true, then replace
%       the background)
%       > 'X': shift the raw data distribution using the log-normal
%       property. Resulting data re LOG-NORMALly distriubted. Default
%       dynamic range of data is [0,255]
%       > 'Y': mean-variance equalize the LOG of the raw data s.t.
%       resulting data re NORMALLY distriubted. Default dynamic range is
%       [-3,3]
%   • 'medfilt': logical indicating whether to median filter the data or
%   not [DEFAULT = 0]
%   • 'labelInfo': structure which is output from labelMIMSstruct.m
%   function containing the 
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author            E-mail                      Version
%  31 July 2015   Amanda Balderrama amanda.gaudreau@gmail.com     0
%   3 Aug  2015   Amanda Balderrama amanda.gaudreau@gmail.com     1
%  28 Sept 2015   Amanda Balderrama amanda.gaudreau@gmail.com     2
%     For tissue foregound segementation, that case where the element
%     approximates an exponential well is considered
%   1 Nov  2015   Amanda Balderrama amanda.gaudreau@gmail.com     2.1
%     Incorporated ability to input hdrtxt and labelMIMSstruct output as
%     well as use the function to define mask

%% ---- Read data from text files, assumes comma delimiter ----
signal_point = 1;
compdir = loaddirfun;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if isstr(data)
  [data,hdrtxt] = readfiles('dat',data,varargin{:});
elseif strmatch('hdrtxt',PropertyNames)
  hdrtxt = PropertyVal{strmatch('hdrtxt',PropertyNames)};
else 
  hdrtxt = [];
end

if strmatch('repback',PropertyNames)
  repback = PropertyVal{strmatch('repback',PropertyNames)};
else
  repback = 1;
end

if strmatch('rotdat',PropertyNames)
  rotdat = PropertyVal{strmatch('rotdat',PropertyNames)};
else
  rotdat = 1;
end

if strmatch('labelstr',PropertyNames)
  labelstr = PropertyVal{strmatch('labelstr',PropertyNames)};
else
  labelstr = [];
end

mimsInfo = MIMSsummary;
di = strmatch(data.dataset,mimsInfo.dsetInfo(:,1));
if rotdat & ~isempty(di)
  rotval = mimsInfo.rotvec(di);
else
  rotval = 0;
end

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 0;
end

if strmatch('saveimgs',PropertyNames)
  saveimgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  saveimgs = 0;
end

if strmatch('outtag',PropertyNames)
  outtag = PropertyVal{strmatch('outtag',PropertyNames)};
else
  outtag = 'rXY';
end

if strmatch('medfilt',PropertyNames)
  mdflt = PropertyVal{strmatch('medfilt',PropertyNames)};
else
  mdflt = 0;
end

if strmatch('labelInfo',PropertyNames)
  labelInfo = PropertyVal{strmatch('labelInfo',PropertyNames)};
else
  labelInfo = [];
end

%>>> TO RUN AS SCRIPT:
% mimsInfo = MIMSsummary;
% dstr = mimsInfo.dsetInfo{38,1};
% rotval = mimsInfo.dsetInfo{38,2};
% [data,hdrtxt] = readfiles('dat',dstr,'labelstr','brainHIPP');
% labelstr = 'brainHIPP';
% useElstr = {'P_','Fe','Cu63','C_13','Zn64','Gd'};

xpth = [0.5,0.0005];
xmapvals = [50,254.5];
sigXtarget =  (log(xmapvals(1)) - log(xmapvals(2)))/(sqrt(2)*(erfinv(2*xpth(1) - 1) - erfinv(1-2*xpth(2))));
muXtarget = log(xmapvals(1)) - sqrt(2)*sigXtarget*erfinv(2*xpth(1)-1); %= log(xmapvals(2)) - sqrt(2)*sigXtarget*erfinv(1-2*xpth(2));

if ~isempty(labelstr)
  [~,hdrtxt,d_div,line_types,type_rows] = divMIMSstruct(data,hdrtxt);
  %[d,hdrtxt,d_div,line_types,type_rows,~,~,~,dchmask] = labelMIMSstruct(d,hdrtxt,'ncc',1,'plotimgs',1);
  c = find(cellfun(@(x) ~isempty(x),strfind(line_types,labelstr)));
  dstr = d_div{c};
else
  c = 1;
  dstr = data;
  type_rows{1} = 1:size(dstr.Time,1);
end

if isempty(labelInfo)
  htlab = hdrtxt([1;type_rows{c}(:)+1],:);
  [~,~,~,~,~,labelInfo] = labelMIMSstruct(dstr,htlab);
  c=1;
end

delonly = rmfield(dstr,skip_fields);
elnames = fieldnames(delonly);
dsize(3) = length(elnames);
[dsize(1),dsize(2)] = size(dstr.Time);
dmat = reshape(struct2array(delonly),dsize(:)');
nanX = sum(isnan(dmat),3); % summed logical image which ranges between 0 - #el indicating how many channels have a nan value for a given pixel
nanXcol = sum(nanX,1); % summed columnwise: 0 - #rows indicating # nan values in a given column
keepcol = nanXcol < dsize(1)*dsize(3)*.3; % keep columns with 70% or more non-nan values
nanXrows = sum(nanX,2); % summed rowwise: 0 - #cols indicating # nan values in a given row
keeprows = nanXrows < dsize(2)*dsize(3)*.3;

tsuMask = labelInfo.tissueMask{c};

dmat = dmat(keeprows,keepcol,:);
tsuMask = tsuMask(keeprows,keepcol,:);
useEli = interp_data_elemrng(dstr,useElstr); %1:dsize(3);%[1     3     4     6     8    10    12    13    14    16    17    18    19];%find(median(dHmat) > 0);%
elnames = elnames(useEli);
duse = dmat(:,:,useEli);

% >>> Perform necessary processing on multichannel image
% --- 1) ROTATING: Rotate image if necessary (because of different
% cassette alignment in ablation block)
Iraw = rotMIMS(duse,rotval);
dsize = size(Iraw);
if length(dsize) < 3
  dsize(3) = 1;
end
tsuMask = logical(rotMIMS(tsuMask,rotval));
Trot = rotMIMS(data.Time,rotval);
maskcc = regionprops(tsuMask,'Image','Orientation','boundingbox',...
  'MajorAxisLength','MinorAxisLength','extent','Eccentricity','centroid');

% --- 2) MEAN VARIANCE EQUALIZE
% --- Flatten image
dflat = reshape(Iraw(:),[prod(dsize(1:2)),dsize(3)]);
dnon0 = dflat > 0;
X = dflat(tsuMask,:);
Xnon0 = X > 0;
Y = log(X+eps);
Y(~Xnon0) = nan;
muYinit = nanmean(Y);
sigYinit = nanstd(Y);
sigYinit(sigYinit > 1) = 1;
CY = zeros(1,dsize(3));
muY = zeros(1,dsize(3));
sigY = zeros(1,dsize(3));
nsig = 4;
Ynsiglims = [(muYinit-sigYinit.*nsig);(muYinit+sigYinit.*nsig)];
% --- Find necessary statistics for mean-variance equalization
for el = 1:dsize(3)
  %disp(el);
  Ypcurve = Y(Y(:,el) >= Ynsiglims(1,el) & Y(:,el) <= Ynsiglims(2,el),el);
  yrng = Ynsiglims(2,el) - Ynsiglims(1,el);
  M = 200;
  bw100 = sigYinit(el)./10;
  yvals = Ynsiglims(1,el):yrng/(M-1):Ynsiglims(2,el);
  [uyvals,ia,ic] = unique(Ypcurve);
  if length(uyvals) < M
    yvals = uyvals';
  end
  hy = hist(Ypcurve,yvals);
  hykde = ksdensity(Ypcurve,yvals,'bandwidth',bw100);
  %figure; plot(yvals,[hy;hykde./max(hykde).*max(hy)]);
  
  yhat = hy;
  niter = 10;
  Xiter = zeros(niter,3);
  A = zeros(3,3);
  C = zeros(niter,1);
  muemp = zeros(niter,1);
  sig2emp = zeros(niter,1);
  
  if length(yvals) < 50
    [~,iymax] = max(hy);
    muemp(niter) = yvals(iymax);
    sig2emp(niter) = mean(abs(diff(yvals)))^2;
    CY(el) = hy(iymax);%C(niter);
    muY(el) = muemp(niter);
    sigY(el) = sqrt(sig2emp(niter));
  else
    for k = 1:niter
      if k == 1
        yiter = yhat;
        useyi = yiter > 1;
      else
        yiter = exp(Xout(1)+Xout(2).*yvals + Xout(3).*yvals.^2);
      end
      if ~any(isinf(yiter))
        y2iter = yiter.^2;
        A(1,:) = [sum(y2iter(useyi)), sum(yvals(useyi).*y2iter(useyi)), sum(yvals(useyi).^2.*y2iter(useyi))];
        A(2,:) = [A(1,2), A(1,3), sum(yvals(useyi).^3.*y2iter(useyi))];
        A(3,:) = [A(2,2), A(2,3), sum(yvals(useyi).^4.*y2iter(useyi))];
        B = [sum(y2iter(useyi).*log(yhat(useyi)));...
          sum(yvals(useyi) .* y2iter(useyi) .* log(yhat(useyi)));...
          sum(yvals(useyi).^2 .* y2iter(useyi) .* log(yhat(useyi)))];
        Xout = linsolve(A,B);
        Xiter(k,:) = Xout(:)';
        C(k) = exp(Xout(1) - Xout(2)^2/(4*Xout(3)));
        muemp(k) = -Xout(2)/(2*Xout(3));
        sig2emp(k) = -1/(2*Xout(3));
      else
        [~,iymax] = max(hy);
        muemp(niter) = yvals(iymax);
        sig2emp(niter) = std(Ypcurve)^2;% mean(abs(diff(hy)))^2;
        CY(el) = C(niter);
        muY(el) = muemp(niter);
        sigY(el) = sqrt(sig2emp(niter));
      end
    end
    
    if muemp(niter) > Ynsiglims(2,el)
      [~,maxhyi] = max(hy);
      muY(el) = yvals(maxhyi);
      CY(el) = hy(maxhyi);
      sigY(el) = sqrt(((yvals(maxhyi+1) - yvals(maxhyi))^2)/(2*log(hy(maxhyi)/hy(maxhyi+1))));
    else
      CY(el) = C(niter);
      muY(el) = muemp(niter);
      if sig2emp(niter) < 0
        sigY(el) = sqrt(0.1);
      else
        sigY(el) = sqrt(sig2emp(niter));
      end
    end
  end
end

if mdflt
  mfsize = [max(round(dsize(1)/150),3),max(round(dsize(2)/150),3)];
else
  mfsize = nan;
end

if ismember('r',outtag)
  doutmat = Iraw;
  if repback
    bkgndval = zeros(1,dsize(3));
    for el = 1:dsize(3)
      bkgndval(el) = Ynsiglims(1,el);
      dflat(~tsuMask,el) = bkgndval(el);
    end
    doutmat = reshape(dflat,dsize);
  else
    bkgndval = nan;
  end
  dout.Iraw = mat2dstruct(doutmat,data,mfsize,elnames);
end

% --- Final MVE code
targetsig = 1.*ones(1,dsize(3));
targetmu = 0.*ones(1,dsize(3));
alphaY = sqrt(targetsig./sigY.^2);
betaY = targetmu - alphaY .* muY;
IYd = repmat(alphaY,[prod(dsize(1:2)),1]).*log(dflat) + repmat(betaY,[prod(dsize(1:2)),1]);
IYd(~dnon0) = -3;
for el = 1:dsize(3)
 IYd(dnon0(:,el) == 0,el) =  min(IYd(dnon0(:,el),el));
end
% --- Replace background with tissue mean (== targetmu)
if repback
  bkgndval = norminv(0.05,targetmu(1),targetsig(1)).*ones(1,dsize(3));
  newbkgnd = repmat(bkgndval,[sum(sum(~tsuMask)),1]);% randn(size(dbkgnd)).*repmat(sigbkgnd,[size(dbkgnd,1),1]) + repmat(mubkgnd,[size(dbkgnd,1),1]);
  IYd(~tsuMask,:) = newbkgnd;
else
  bkgndval = nan;
end
doutmat = reshape(IYd,dsize);

if plotimgs
  IYmve = doutmat;
end

if ismember('Y',outtag)
  dout.IYmve = mat2dstruct(doutmat,data,mfsize,elnames);
end


% --- Final X shifted according to xmapvals parameter (beginning of script)
targetsig = sigXtarget.*ones(1,dsize(3));
targetmu = muXtarget.*ones(1,dsize(3));
alpha = sqrt(targetsig.^2./sigY.^2);%sigYhat);%
beta = targetmu - alpha .* muY;%muYhat;%
IXshd = repmat(alpha,[prod(dsize(1:2)),1]).*log(dflat) + repmat(beta,[prod(dsize(1:2)),1]);
% --- Replace background with tissue mean (== targetmu)
if repback
  bkgndval = norminv(0.05,targetmu(1),targetsig(1)).*ones(1,dsize(3));
  newbkgnd = repmat(bkgndval,[sum(sum(~tsuMask)),1]);% randn(size(dbkgnd)).*repmat(sigbkgnd,[size(dbkgnd,1),1]) + repmat(mubkgnd,[size(dbkgnd,1),1]);
  IXshd(~tsuMask,:) = newbkgnd;
else
  bkgndval = nan;
end
IXshd = exp(IXshd);
%IXshd(~dnon0) = 0;
doutmat = reshape(IXshd,dsize);

if plotimgs
  IXshift = doutmat;
end

if ismember('X',outtag)
  dout.IXshift = mat2dstruct(doutmat,data,mfsize,elnames);
end

mveparams.dataset = data.dataset;
mveparams.elnames = elnames;
mveparams.tsuMask = tsuMask;
mveparams.dsize = dsize;
mveparams.outimgstr = outtag;
mveparams.labelstr = labelstr;
mveparams.bkgndval = bkgndval;
mveparams.rotateval = rotval;
mveparams.medfilter = mfsize;
mveparams.Ymu0 = muYinit;
mveparams.Ysig0 = sigYinit;
mveparams.Yrng = Ynsiglims;
mveparams.Ymucalc = muY;
mveparams.Ysigcalc = sigY;
dout.mveparams = mveparams;

if plotimgs
  %d3allLab = zeros(prod(dsize(1:2)),3);
  %or ch = 1:size(labelInfo.dchMask,3)
  %  d3allLab(tsuMask,:) = repmat([1,0,1],tsuMask(:),1);
  %end
  %d3allLab = reshape(d3allLab,[dsize(1),dsize(2),3]);
  d3allLab = intensity2rgb(abs(tsuMask-1),'c','i2nan',0);%
  sp1 = [0.1300    0.5838    0.3347    0.3412];
  sp2 = [0.5703    0.5838    0.3347    0.3412];
  sp3 = [0.1300    0.1100    0.3347    0.3412];
  sp4 = [0.5703    0.1100    0.3347    0.3412];
  bkgdlevel = 0.15;
  fsz = 20;

  for el = 1:dsize(3)
    Ypcurve = Y(Y(:,el) >= Ynsiglims(1,el) & Y(:,el) <= Ynsiglims(2,el),el);
    yrng = Ynsiglims(2,el) - Ynsiglims(1,el);
    M = 200;
    bw100 = sigYinit(el)./10;
    yvals = Ynsiglims(1,el):yrng/(M-1):Ynsiglims(2,el);
    [uyvals,ia,ic] = unique(Ypcurve);
    if length(uyvals) < M
      yvals = uyvals';
    end
    modelyvals = yvals(1):yrng/(M-1):yvals(end);
    hy = hist(Ypcurve,yvals);
    hykde = ksdensity(Ypcurve,yvals,'bandwidth',bw100);
    
    yhat = hy;
    useyi = yhat > 0;
    nbins = 50;
    figure('name',sprintf('%slnnormfit_%s',labelstr,elnames{el}),...sprintf('%s %s MVE comparison',evalf,elnames{f,useEl(el)}),...
      'position',[9   339   900   657]);%[49         155        1008         795]);%
    %---------------------------------------------
    subplot(2,2,1); im2 = imagesc(Iraw(:,:,el),exp([yvals(1), yvals(end)]));colormap(gray);hold on;% auto_thresh(X(:,el))); colormap(gray);hold on; 
    cbh = colorbar;%('location','southoutside');
    im1 = image(d3allLab); 
    set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*bkgdlevel);
    title('Raw [cps] = X','fontsize',fsz);%,'interpreter','latex');
    set(gca,'position',sp1,'xticklabel','','yticklabel','');
    set(cbh,'position',[sp1(1)+sp1(3),sp1(2)+.015,0.01,sp1(4)-.03],'fontsize',fsz)
    %---------------------------------------------
    logdimg = log(Iraw(:,:,el)+eps); logdth = auto_thresh(logdimg);
    subplot(2,2,2); imagesc(logdimg,[yvals(1), yvals(end)]); colormap(gray); hold on;
    cbh = colorbar;%('location','southoutside');
    im1 = image(d3allLab); 
    set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*bkgdlevel);
    set(gca,'position',sp2,'xticklabel','','yticklabel','');
    set(cbh,'position',[sp2(1)+sp2(3),sp2(2)+.015,0.01,sp2(4)-.03],'fontsize',fsz)
    title('ln(Raw [cps]) = ln X = Y','fontsize',fsz);%,'interpreter','latex');
    %---------------------------------------------
    subplot(2,2,4);
    plot(yvals(useyi),yhat(useyi),'o')
    hold on; plot(yvals,hykde./max(hykde).*max(yhat),'c');
    yfixedC = CY(el).*exp(-(yvals(:) - muY(el)).^2./(2*sigY(el)^2));% normpdf(yvals,muemp(niter),sqrt(sig2emp(niter))).*C(niter);
    plot(yvals,yfixedC,'m');
    title({[sprintf('%s Tissue Values, ln[cps]',strrep(elnames{el},'_','\_'))],...
      [sprintf('\\mu = %1.2f, \\sigma = %1.3f',...
      muY(el),sigY(el))]},'fontsize',fsz)
    legend({'Raw Hist','KDE Fit','Normal Fit'})
    xlabel('Y = ln X','fontsize',fsz);%,'interpreter','latex')
    ylabel('Frequency of occurrence','fontsize',fsz)
    axis tight; grid on;
    set(gca,'position',sp4,'fontsize',fsz);
    %---------------------------------------------
    [hcnts,xcnts] = hist(X(:,el),exp(yvals));
    subplot(2,2,3); plot(xcnts,[hcnts]','o'); hold on;
    plot(xcnts,(yfixedC./max(yfixedC)).*max(hcnts),'m');
    %hN = plot(xcnts,(yhat./max(yhat)).*max(hcnts),'g');
    title({[sprintf('%s Tissue Values',strrep(elnames{el},'_','\_'))],...
      [sprintf('')]},'fontsize',fsz)
    legend({'Raw Hist','Log-Norm Fit'})
    xlabel('X','fontsize',fsz);%,'interpreter','latex')
    ylabel('Frequency of occurrence','fontsize',fsz)
    xlim([exp(yvals(1)),exp(yvals(end-1))]); grid on;
    set(gca,'position',sp3,'fontsize',fsz);
    
    figure('name',sprintf('%sshift_%s',labelstr,elnames{el}),...sprintf('%s %s MVE comparison',evalf,elnames{f,useEl(el)}),...
      'position',[9   369   900   657]);%[9   339   963   657]);%[49         155        1008         795]);
    %---------------------------------------------
    subplot(2,2,1); im2 = imagesc(IXshift(:,:,el),[0,255]); colormap(gray); hold on;
    %subplot(2,2,1); image(intensity2rgb(IXshift(:,:,el),'w','thresh',[0,255]));
    cbh = colorbar;%('location','southoutside');
    im1 = image(d3allLab); hold on;
    set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*bkgdlevel);
    title(sprintf('X_{shifted} ~ ln [N(%1.2f,%1.3f)]',muXtarget,sigXtarget),'fontsize',fsz);%,'interpreter','latex');
    set(gca,'position',sp1,'xticklabel','','yticklabel','');
    set(cbh,'position',[sp1(1)+sp1(3),sp1(2)+.015,0.01,sp1(4)-.03],'fontsize',fsz)
    %---------------------------------------------
    subplot(2,2,2); im2 = imagesc(IYmve(:,:,el),[-3,3]); colormap(gray); hold on;
    cbhb = colorbar;%('location','southoutside');
    im1 = image(d3allLab); hold on;
    set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*bkgdlevel);
    title(sprintf('Y_{MVE} ~ N(0,1)'),'fontsize',fsz);%,'interpreter','latex');
    set(gca,'position',sp2,'xticklabel','','yticklabel','');
    set(cbhb,'position',[sp2(1)+sp2(3),sp2(2)+.015,0.01,sp1(4)-.03],'fontsize',fsz)
    %---------------------------------------------
    dmin = -3;%min(IYd(Hmaskmat(:,:,4),el));
    dmax = 0;%max(IYd(Hmaskmat(:,:,4),el));
    uxvals = 0;
    while length(uxvals) < 10
      dmax = dmax + 3;
      xrng = dmax - dmin;
      xvals = dmin:(xrng)/(M-1):dmax;
      IYduse = IYd(tsuMask & ...
        reshape(IYd(:,el) >= dmin,dsize(1),dsize(2)) & ...
        reshape(IYd(:,el) <= dmax,dsize(1),dsize(2)),el);
      [uxvals,ia,ic] = unique(IYduse);
    end
    
    if length(uxvals) < M
      xvals = uxvals';
    end
    modelxvals = xvals(1):xrng/(M-1):xvals(end);
    [hcnts,xcnts] = hist(IYd(tsuMask,el),xvals);
    subplot(2,2,4); plot(xcnts(hcnts > 0),[hcnts(hcnts > 0)]','o'); hold on;
    cp = get(gca,'position');
    yN = normpdf(modelxvals);% mean(lg10dtsu(usei,el)),std(lg10dtsu(usei,el)));
    if modelxvals(1) ~= 0
      hN = plot(modelxvals,(yN./max(yN).*max(hcnts(2:(end-1)))),'m');
      axis tight; 
      ylim([min(hcnts),max([min(hcnts) + 1,hcnts(2:end-1)])])
    else
      hN = plot(modelxvals,(yN./max(yN).*max(hcnts(1:(end-1)))),'m');
      axis tight; 
      ylim([min(hcnts),max([min(hcnts) + 1,hcnts(1:end-1)])])
    end
    %yhat = normpdf(xN,muYhat(el),sigYhat(el));
    %plot(xN,(yhat./max(yhat)).*max(hcnts),'g');
    title(sprintf('%s MVE Tissue Values',strrep(elnames{el},'_','\_')),'fontsize',fsz)
    xlabel('Y_{MVE}','fontsize',fsz);%,'interpreter','latex')
    ylabel('Frequency of occurrence','fontsize',fsz)
    set(gca,'position',sp4,'fontsize',fsz);
    %---------------------------------------------
    dmin = log(1);%log(min(IXshd(Hmaskmat(:,:,4),el)));
    dmax = 0;
    uxvals = 0;
    while length(uxvals) < 10
      dmax = dmax + log(xmapvals(2));
      xrng = dmax - dmin;
      xvals = dmin:(xrng)/(M-1):dmax;
      IXshduse = IXshd(tsuMask & ...
        reshape(IXshd(:,el) > 0,dsize(1),dsize(2)) & ...
        reshape(IXshd(:,el) <= exp(dmax),dsize(1),dsize(2)),el);
      [uxvals,ia,ic] = unique(IXshduse);
    end
    
    if length(uxvals) < M
      xvals = [0,log(uxvals')];
    end
    modelxvals = xvals(1):xrng/(M-1):xvals(end);
    yN = normpdf(modelxvals,muXtarget,sigXtarget);% mean(lg10dtsu(usei,el)),std(lg10dtsu(usei,el)));
    [hcnts,xcnts] = hist(IXshd(tsuMask,el),exp(xvals));
    subplot(2,2,3); plot(xcnts(hcnts > 0),[hcnts(hcnts > 0)]','o'); hold on;
    hN = plot(exp(modelxvals),(yN./max(yN).*max(hcnts(2:(end-1)))),'m');
    title(sprintf('%s Shifted Tissue Values',strrep(elnames{el},'_','\_')),'fontsize',fsz)
    xlabel('X_{shifted}','fontsize',fsz);%,'interpreter','latex')
    ylabel('Frequency of occurrence','fontsize',fsz)
    axis tight;hold off;ylim([0,min(max(hcnts),max(hcnts(2:(end-1)))*1.5)])
    set(gca,'position',sp3,'fontsize',fsz);
  end
end

if saveimgs
  save_open_figures(strcat(compdir.lafldr,data.dataset,'\Images\'),[],[],...
    strcat(data.dataset));
end

%=========================================================================
% SUBFUNCTIONS
%=========================================================================

  function ds = mat2dstruct(doutmat,data,mfsize,elnames)
    ds = struct();
    if ~isnan(mfsize)
      for k = 1:length(elnames)
        ds = setfield(ds,elnames{k},medfilt2(doutmat(:,:,k),mfsize));
      end
    else
      for k = 1:length(elnames)
        ds = setfield(ds,elnames{k},doutmat(:,:,k));
      end
    end
    ds.Time = Trot;
    ds.line_datevec = data.line_datevec;
    ds.dataset = data.dataset;
    fldorder = cell(1,length(elnames) + length(skip_fields));
    fldorder(1:length(elnames)) = sort(elnames);
    fldorder((length(elnames)+1):end) = skip_fields;
    ds = orderfields(ds,fldorder);
  end
end