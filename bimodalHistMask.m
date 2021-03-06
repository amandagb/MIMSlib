function [Imask,varargout] = bimodalHistMask(I,varargin)
%% function: calibrationMask
% Description:
% Example:
% c = find(ismember(line_types,'brainHIPP'));
% tsuel = 'S_34';
% dmaskraw = getfield(d_div{c},tsuel);
% [~,~,dthmask] = auto_thresh(dmaskraw);
% % obj = gmdistribution.fit(dthmask(:),2);
% tissueMask = bimodalHistMask(dthmask,'Nmedfilter',3,'sqmask',0,'quantfrac',0.01,...
%   'nbins',50,'plotimgs',1,'fill',0,'close',1,'thupperI',0,'ncc',1);
%
% c = find(ismember(line_types,'Calib_5'));
% laoffel = 'C_13';
% dmaskraw = getfield(d_div{c},laoffel);
% [~,~,dthmask] = auto_thresh(dmaskraw);
% laoff = bimodalHistMask(dthmask,'Nmedfilter',1,'sqmask',0,'quantfrac',0.01,...
%   'nbins',50,'plotimgs',1,'fill',1,'close',1,'thupperI',0,'ncc',1);
%
% TO REPLICATE v1 behavior: set ----   'usenorm', 0
% INPUTS ----------------------------------------------------------------
% I:    Input image (M x N)
% varargin - 'PropertyName','PropertyValue'
%   � 'Imasks':   structure (see outputs) => if indicates, F/B segmentation
%   procedure will not be executed and the masks in the structure will be
%   used instead
%   � 'Nmedfilter': integer (best if odd) indicating the width of the median
%   filter used to smooth the data [DEFAULT = 9]
%   � 'nbins':  integer indicating how many histogram bins to use for
%   binning the data. IF nbins == 0, then unique values will be used
%   instead of a number of histogram bins [DEFAULT = 50]
%   � 'usenorm': logical indicating whether to use gmm algorithm or not
%   [DEFAULT = 1]
%   � 'minpnt': logical indicating whether threshold the data such that the
%   the value corresponds to the point with minimum bin occupancy (minimum
%   point in histogram) [DEFAULT = 0]
%   � 'Ngmm': number of components in gmm [DEFAULT = 2]
%   � 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   � 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   � 'dataname': string indicating the data's label
%   � 'sqmask': logical indicating whether to square the mask or not
%   [DEFAULT = 0]
%   � 'fill': logical indicating whether to fill the mask or not [DEFAULT = 1]
%   � 'ncc': integer which indicates the number of connected components to
%   be enforced in the image [DEFAULT = [] --- no number enforced]
%
% OUTPUTS ---------------------------------------------------------------
% Imask:
% varargout
%   {1} structOut: Output structure with fields containing data used to
%   compute peak locations and threshold
%
%
%  Date           Author              E-mail                      Version
%   1 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
%   3 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%     Estimate the GMM parameters & cut off the foreground intensity based
%     on
%  24 June 2015   Amanda Balderrama   amandagbalderrama@gmail.com     3
%  28 Sept 2015   Amanda Balderrama   amanda.gaudreau@gmail.com       4
%     For tissue foregound segementation, that case where the element
%     approximates an exponential well is considered

compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 0;
end

if strmatch('axh',PropertyNames)
  ah = PropertyVal{strmatch('axh',PropertyNames)};
else
  ah = [];
end

if strmatch('figh',PropertyNames)
  fh = PropertyVal{strmatch('figh',PropertyNames)};
else
  fh = [];
end

if strmatch('saveimgs',PropertyNames)
  saveimgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  saveimgs = 0;
end

if strmatch('sqmask',PropertyNames)
  sqMask = PropertyVal{strmatch('sqmask',PropertyNames)};
else
  sqMask = 0;
end

if strmatch('quantfrac',PropertyNames) % quantile fraction
  qth = PropertyVal{strmatch('quantfrac',PropertyNames)};
else
  qth = 0.05;
end

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 50;
end

if strmatch('usenorm',PropertyNames)
  usenorm = PropertyVal{strmatch('usenorm',PropertyNames)};
else
  usenorm = 1;
end

if strmatch('minpnt',PropertyNames)
  minpnt = PropertyVal{strmatch('minpnt',PropertyNames)};
else
  minpnt = 0;
end

if strmatch('ngmm',PropertyNames)
  k = PropertyVal{strmatch('ngmm',PropertyNames)};
else
  k = 2;
end

if strmatch('fill',PropertyNames)
  maskfill = PropertyVal{strmatch('fill',PropertyNames)};
else
  maskfill = 1;
end

if strmatch('close',PropertyNames)
  maskclose = PropertyVal{strmatch('close',PropertyNames)};
else
  maskclose = 0;
end

if strmatch('ncc',PropertyNames)
  numcc = PropertyVal{strmatch('ncc',PropertyNames)};
else
  numcc = [];
end

if strmatch('nmodes',PropertyNames)
  npks = PropertyVal{strmatch('nmodes',PropertyNames)};
else
  npks = 2;
end

if strmatch('thupperI',PropertyNames)
  thupperI = PropertyVal{strmatch('thupperI',PropertyNames)};
else
  thupperI = 0;
end

if strmatch('dataname',PropertyNames)
  fname = PropertyVal{strmatch('dataname',PropertyNames)};
else
  fname = '';
end

ubins = unique(I(:));
if nbins && length(ubins) > nbins
  [hcnts,xcnts] = hist(I(:),nbins);
else
  xcnts = ubins;
  hcnts = hist(I(:),xcnts);
end

if strmatch('nmed',PropertyNames)
  nMF = PropertyVal{strmatch('nmed',PropertyNames)};
else
  nMF = max([length(hcnts)/20,3]);
end

MFUcnts = medfilt1(hcnts,nMF);
[M,N] = size(I);

%-----------------------------------------
% Visualization of histogram
if plotimgs
  if isempty(fh)
    fh = figure('name',fname);
  end
  plot(xcnts,[hcnts;MFUcnts]'); hold on; axis tight;
  xlabel('counts');
  ylabel('# occurrences in image');
  title(sprintf('%s pixel count distribution',strrep(fname,'_','\_')))
end
%-----------------------------------------

%-----------------------------------------
% Determine whether there are peaks in the histogram -- detect whether it's
% more of an exponential decay shape
%-----------------------------------------

if length(MFUcnts) == npks
  ppindord = [1:npks];
  pkpos = [1:npks];
  pkval = MFUcnts(1:npks);
  upperIco = xcnts(end);
  lowerIco = mean(xcnts);
else
  [pks,pind] = findpeaks(MFUcnts);
  if isempty(pks) % this occurs if the only peak is the peak at the beginning of the histogram
    pks = [MFUcnts(1),MFUcnts(end)];
    pind = [1,length(MFUcnts)];
  else
    pks = [MFUcnts(1),pks(:)'];
    pind = [1,pind(:)'];
    %pind = pind - 1;
  end
  
  [pksord,indord] = sort(pks,'descend');
  if length(pksord) > (npks + 1)
    % This code block identifies regions to "clump" together
    keepi = 1;
    nvalid = 0;
    dpks = abs(diff(pind(indord)));%diff(pksord);
    pntco = 5;
    for di = 1:(length(pind)-1)
      ph = pksord(di:di+1);
      piends = [pind(indord(di+1)), pind(indord(di))];
      dpi = abs(diff(piends));
      if abs(dpi) < pntco && nvalid < npks
        btwi = min(piends):max(piends);
        mv = min(MFUcnts(btwi));
        if mv < min(ph)*0.5
          keepi = [keepi,di+1];
          nvalid = nvalid + 1;
        elseif ismember(di,keepi)
          pksord(di+1) = pksord(di);
          indord(di+1) = indord(di);
        end
      else
        keepi = [keepi,di+1];
        nvalid = nvalid + 1;
      end
    end
    dpks_large = keepi;
  elseif MFUcnts(1) < min(pks) && length(pksord) > npks
    dpks_large = 2:length(pksord);%1: npks;%
  else
    dpks_large = 1:npks;
  end
  
  if length(dpks_large) < npks
    dpks_large = 1:npks;
  end
  ipnts = [0,dpks_large(1:npks)];
  % pkpos = [pind(indord(1:ipnts(1))),pind(indord(ipnts(1)+1:ipnts(2)))
  pkpos = zeros(1,npks);
  pkval = zeros(1,npks);
  if length(pind) == npks
    for i = 1:npks
      pkpos(i) = pind(  indord(  ipnts(i+1 )  )  );
      pkval(i) = pks(indord(ipnts(i+1)));
      %pkpos(i) = round(mean(pind(indord((ipnts(i)+1):ipnts(i+1)))));
      %pkval(i) = mean(pks(indord((ipnts(i)+1):ipnts(i+1))));
    end
  else
    for i = 1:npks
      pkpos(i) = pind(  indord(  ipnts(  min([i+1,length(pind)])  )  )  );
      pkval(i) = pks(indord(ipnts(min([i+1,length(pind)]))));
      %pkpos(i) = round(mean(pind(indord((ipnts(i)+1):ipnts(i+1)))));
      %pkval(i) = mean(pks(indord((ipnts(i)+1):ipnts(i+1))));
    end
  end
  
  [ppord,ppindord] = sort(pkpos,'ascend');
  if any(diff(ppord) <= 1 )
    closeind = [1,find(diff(ppord) ~= 1)+1];
    pkpos = pkpos(closeind);
    [ppord,ppindord] = sort(pkpos,'ascend');
  end
  
  pkpos1 = pkpos(ppindord(1));
  pkval1 = pkval(ppindord(1));
  if length(ppindord) < 2
    lowerIco = pkval1;
    upperIco = pkval1 + abs(diff(xcnts(1:2)));
  else
    pkpos2 = pkpos(ppindord(2));
    pkval2 = pkval(ppindord(2));
    Ival1 = xcnts(pkpos1);
    Ival2 = xcnts(pkpos2);
    if plotimgs
      ylim([0,pkval2*1.5])
      plot(xcnts([pkpos1;pkpos1]),[0,pkval1],':k')
      plot(xcnts([pkpos2;pkpos2]),[0,pkval2],':k')
    end
    
    if qth && ~usenorm
      cntco = 0;
    else
      [cntco,bestind] = min(MFUcnts(pkpos1:pkpos2));
    end
    TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
    medbtwnval = median(MFUcnts(pkpos1:pkpos2));
    %     ccTF = bwconncomp(TFvec);
    %     while ccTF.NumObjects > 1
    %       cntco = cntco + 1;
    %       bestind = find(MFUcnts(pkpos1:pkpos2) == cntco);
    %       bestind = bestind(1);
    %       TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
    %       ccTF = bwconncomp(TFvec);
    %     end
    
    if plotimgs
      hco = plot([xcnts(1),xcnts(end)],[cntco,cntco],'--r');
      htf = plot(xcnts(pkpos1:pkpos2),abs(TFvec-1).*pkval2,'g');
    end
    
    if minpnt
      lowerIco = xcnts(bestind+pkpos1-1);
    elseif usenorm
      if qth <= 0
        qth = 0.01;
      end
      ccTF.NumObjects = 5;
      qth = qth - 0.01;
      while ccTF.NumObjects > 1 && qth < 0.25 && qth >= 0
        qth = qth + 0.01;
        %   Igmm = I;
        %   Igmm(I <= xcnts(bestind)) = 1;
        %   Igmm(I > xcnts(bestind)) = 2;
        %   obj = gmdistribution.fit(I(isnan(Igmm) == 0),2,'start',Igmm(isnan(Igmm) == 0));
        %   if plotimgs
        %     pxvals = pdf(obj,xcnts(:));
        %     plot(xcnts(:),pxvals./max(pxvals).*max(hcnts),'c');
        %   end
        Itsu = I(I > xcnts(bestind+pkpos1-1));
        mu = mean(Itsu(:));
        sig = std(Itsu(:));
        lowerIco = max(norminv(qth,mu,sig),xcnts(bestind+pkpos1-1));
        [val,coind] = min(abs(xcnts - lowerIco));
        cntco = MFUcnts(coind);
        TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
        upperIco = norminv(1-qth,mu,sig);
        if plotimgs
          delete(hco)
          delete(htf)
          if exist('hN')
            delete(hN)
          end
          %hco = plot([0;length(xcnts)],[cntco,cntco],'--r');
          hco = plot(xcnts([1,end]),[cntco,cntco],'--r');
          htf = plot(xcnts(pkpos1:pkpos2),abs(TFvec-1).*pkval2,'g');
          Isp = round((upperIco - lowerIco)/100);
          xN = lowerIco:Isp:upperIco;
          yN = normpdf(xN,mu,sig);
          hN = plot(xN,(yN./max(yN)).*pkval2,'m');
        end
        ccTF = bwconncomp(TFvec);
        if cntco > medbtwnval*5 && medbtwnval > 0 || sum(TFvec) == length(TFvec)
          ccTF.NumObjects = 5;
          qth = qth - 0.011;
        end
      end
    else
      if qth
        ccTF.NumObjects = 5;
        while ccTF.NumObjects > 1 && qth < 0.25
          qth = qth + 0.01;
          cntco = quantile(MFUcnts(pkpos1:pkpos2),qth);
          TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
          if plotimgs
            delete(hco)
            delete(htf)
            %hco = plot([0;length(xcnts)],[cntco,cntco],'--r');
            hco = plot(xcnts([1,end]),[cntco,cntco],'--r');
            htf = plot(xcnts(pkpos1:pkpos2),abs(TFvec-1).*pkval2,'g');
          end
          ccTF = bwconncomp(TFvec);
        end
        if plotimgs
          if saveimgs
            export_fig(strcat(cdir.mtngpath,fname,'_photodist'),'-png');
          end
          close(fh);
        end
      end
      TFveci = find(MFUcnts(pkpos1:pkpos2) <= cntco);
      bestind = TFveci(1)+pkpos1; %bestind + pkpos1;%TFveci(end)+pkpos1;% TFveci(round(length(TFveci)/2))+pkpos1; %
      lowerIco = mean(xcnts(TFveci+pkpos1-1));% mean([xcnts(bestind),xcnts(pkpos1)]);%UI(bestind);
      
      TFveci = find(MFUcnts(pkpos2+1:end) <= 2);
      if isempty(TFveci)
        TFveci = length(MFUcnts);
      end
      bestind = TFveci(1)+(pkpos2+1);
      upperIco = xcnts(end);% xcnts(bestind); %---> Upper and lower intensity cut offs are found and applied to the image
    end
  end
end

if ~thupperI
  upperIco = max(I(:));
end

structOut.histcnts = hcnts;
structOut.histcntsMF = MFUcnts;
structOut.histbins = xcnts;
structOut.upperCO = upperIco;
structOut.lowerCO = lowerIco;
structOut.peakI = xcnts(pkpos(ppindord));
structOut.peakhistval = pkval(ppindord);

%====================================================
% Uses a variety of morphological operations to create final tissue
% mask(fill, then open)
%====================================================
maskev = zeros(size(I));
Pmask0 = I > lowerIco & I <= upperIco;
if ~isempty(numcc)
  ccinfo = bwconncomp(Pmask0);
  lencc = cellfun(@(x) length(x),ccinfo.PixelIdxList);
  [valLen,ccind] = max(lencc);
  Pmask0 = false(size(Pmask0));
  Pmask0(ccinfo.PixelIdxList{ccind}) = true;
end

if maskfill
  Pmaskfull = imfill(Pmask0,'holes');
else
  Pmaskfull = Pmask0;
end

se = strel('disk',2,4);
if maskclose
  Pmaskfull = imclose(Pmaskfull,se);
end
% figure; imagesc(I); hold on; contour(Pmaskfull,[0,0],'w')
% figure; imagesc(Pmaskfull);
% figure; plot(sum(Pmaskfull))
if sqMask
  rowsum = sum(Pmaskfull');
  colsum = sum(Pmaskfull);
  maskrows = find(rowsum >= median(rowsum)*.9);
  maskcol = find(colsum >= length(maskrows)*.9);
  Imask = false(size(Pmaskfull));
  Imask(maskrows,maskcol) = Pmaskfull(maskrows,maskcol);
else
  Imask = Pmaskfull;
end
% figure; imagesc(Imask)

varargout{1} = structOut;
end
