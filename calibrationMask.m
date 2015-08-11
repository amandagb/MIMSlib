function [Imask] = calibrationMask(I,varargin)
%% function: calibrationMask
% Description:
% Example:
%     ALLmask = calibrationMask(d4mask,'Nmedfilter',1,'quantfrac',0,'nbins',50,'plotimgs',1,'sqmask',0);
%
% INPUTS ----------------------------------------------------------------
% I:    Input image (M x N)
% varargin - 'PropertyName','PropertyValue'
%   • 'Imasks':   structure (see outputs) => if indicates, F/B segmentation
%   procedure will not be executed and the masks in the structure will be
%   used instead
%   • 'Nmedfilter': integer (best if odd) indicating the width of the median
%   filter used to smooth the data [DEFAULT = 9]
%   • 'nbins':  integer indicating how many histogram bins to use for
%   binning the data. IF nbins == 0, then unique values will be used
%   instead of a number of histogram bins [DEFAULT = 50]
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'dataname': string indicating the data's label
%   • 'sqmask': logical indicating whether to square the mask or not
%   [DEFAULT = 1]
%
% OUTPUTS ---------------------------------------------------------------
% Imask: 
%
%  Date           Author              E-mail                      Version
%  23 Mar  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1

compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

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

if strmatch('sqmask',PropertyNames)
  sqMask = PropertyVal{strmatch('sqmask',PropertyNames)};
else
  sqMask = 1;
end

if strmatch('nmed',PropertyNames)
  nMF = PropertyVal{strmatch('nmed',PropertyNames)};
else
  nMF = 30;
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

if strmatch('dataname',PropertyNames)
  fname = PropertyVal{strmatch('dataname',PropertyNames)};
else
  fname = '';
end

if nbins
  [hcnts,xcnts] = hist(I(:),nbins);
else
  xcnts = unique(I(:));
  hcnts = hist(I(:),xcnts);
end
MFUcnts = medfilt1(hcnts,nMF);
[M,N] = size(I);

%-----------------------------------------
% Visualization of histogram
if plotimgs
  fth = figure('name',fname);
  plot(xcnts,[hcnts;MFUcnts]'); hold on; axis tight;
  xlabel('counts');
  ylabel('# occurrences in image');
  title(sprintf('%s pixel count distribution',strrep(fname,'_','\_')))
end
%-----------------------------------------

[pks,pind] = findpeaks(MFUcnts);
pks = [MFUcnts(1),pks(:)'];
pind = [1,pind(:)'];
[pksord,indord] = sort(pks,'descend');
dpks = diff(pind(indord));%diff(pksord);
dpks_large = find(abs(dpks) > 5);
ipnts = dpks_large(1:2);
% pkpos = [pind(indord(1:ipnts(1))),pind(indord(ipnts(1)+1:ipnts(2)))
pkpos1 = round(mean(pind(indord(1:ipnts(1)))));
pkval1 = mean(pks(indord(1:ipnts(1))));
pkpos2 = round(mean(pind(indord(ipnts(1)+1:ipnts(2)))));
pkval2 = mean(pks(indord(ipnts(1)+1:ipnts(2))));
if pkpos2 < pkpos1
  pp1h = pkpos1;
  pv1h = pkval1;
  pp2h = pkpos2;
  pv2h = pkval2;
  pkpos1 = pp2h;
  pkpos2 = pp1h;
  pkval1 = pv2h;
  pkval2 = pv1h;
end
Ival1 = xcnts(pkpos1);
Ival2 = xcnts(pkpos2);
if plotimgs
  ylim([0,pkval2*1.5])
  plot(xcnts([pkpos1;pkpos1]),[0,pkval1],':k')
  plot(xcnts([pkpos2;pkpos2]),[0,pkval2],':k')
end

if qth
  cntco = 0;
else
  [cntco,bestind] = min(MFUcnts(pkpos1:pkpos2));
end
TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
if plotimgs
  hco = plot([xcnts(1),xcnts(end)],[cntco,cntco],'--r');
  htf = plot(xcnts(pkpos1:pkpos2),abs(TFvec-1).*pkval2,'g');
end

if qth
  ccTF.NumObjects = 5;
  while ccTF.NumObjects > 1
    qth = qth + 0.01;
    cntco = quantile(MFUcnts(pkpos1:pkpos2),qth);
    TFvec = MFUcnts(pkpos1:pkpos2) <= cntco;
    if plotgray
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
    close(fth);
  end
end

TFveci = find(MFUcnts(pkpos1:pkpos2) <= cntco);


bestind = TFveci(1)+pkpos1; %bestind + pkpos1;%TFveci(end)+pkpos1;% TFveci(round(length(TFveci)/2))+pkpos1; %
lowerIco = mean([xcnts(bestind),xcnts(pkpos1)]);%UI(bestind);

TFveci = find(MFUcnts(pkpos2+1:end) <= 2);
if isempty(TFveci)
  TFveci = length(MFUcnts);
end
bestind = TFveci(1)+(pkpos2+1);
upperIco = xcnts(end);% xcnts(bestind); %---> Upper and lower intensity cut offs are found and applied to the image

%====================================================
% Uses a variety of morphological operations to create final tissue
% mask(fill, then open)
%====================================================
maskev = zeros(size(I));
Pmask0 = I > lowerIco & I < upperIco;
Pmaskfull = imfill(Pmask0,'holes');
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
end
