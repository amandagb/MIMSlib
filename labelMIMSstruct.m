function [d,hdrtxt,d_div,line_types,type_rows,outStruct,varargout] = ...
  labelMIMSstruct(d,hdrtxt,varargin)
%% function: labelMIMSstruct
% Description:
% Example:
%     [d,hdrtxt,d_div,line_types,type_rows,labelInfo] = labelMIMSstruct(d,hdrtxt);
%
% INPUTS ----------------------------------------------------------------
% d:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% hdrtxt:
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
%   • 'runVersion': numeric indciating version of function to run. V3 was
%   developed for LacSac experiments, but doesn't work as well as V2
%   developed for MIMS brain experiments (Mar - Sept 2015)
%   • 'tsuFGch': string or cell of strings that can be interpreted by the
%   interp_data_elemrng function which specify the isotopes the user would
%   like to use for masking the tissue from the the rest of the image. IF
%   user would like to specify different "channels" or isotopes for
%   different tissue line types, they should create an nt x 1 cell (where
%   nt is the number of types) of cellstr's which has corresponding strings
%   associated with the line_types mask
%   IF NOT INDICATED, the default method to detect good candidate channels
%   will be used
%   • 'tsula0ch': See 'tsuFGch' descirption. These will instead be used to
%   separate laser off period from laser on period.
%   • 'stndFGch': See 'tsuFGch' description. These will instead be used to
%   separate the "object data" portion of the image corresponding to the
%   laser on period from the laser off period
%   • 'blankout': logical indicating whether to output an arbitrary or
%   blank output versus a true masking function output [DEFAULT = 0]
%   • 'minla0period': numeric indicating the time (in seconds) of the laser
%   off period [DEFAULT = 5]. If there is an inconsistent or unknown laser
%   off period, functionality is best if 'tsula0ch' and/or 'stndFGch' are
%   specified.
%
% OUTPUTS ---------------------------------------------------------------
% d: same as input
% hdrtxt: same as input
% d_div:
% line_types:
% type_rows:
% outStruct: structure with the following fields
%   • labelblocks: ntypes x 1 cell where each cell contains an
%       ntyperows x ntypecols matrix with 5 possible labels
%         {0 = no data; 1 = laser off beginning; 2 = laser off end; 3 =
%         background or imaged area; 4 = tissue}
%   • typerowhdrs = {'Block start row','Block end row','# rows in block','Block type label'};
%   • typerowInfo:
%   • dLabelMask:
%   • laONmask:
%   • tissueMask:
%   • dchMask:
%
%  Date           Author              E-mail                      Version
%   1 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
%  28 Sept 2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%   9 Oct  2015   Amanda Balderrama   amandagbalderrama@gmail.com     3
%     New method for defining masks: instead of voting, a mean log image is
%     formed (which is more robust to outliers) then peaks are detected in
%     the log-cps domain
%  13 Oct  2015   Amanda Balderrama   amandagbalderrama@gmail.com     3.1
%     Merging of new and old methods -- mean log method doesn't work well
%     for MIMS brain datasets
%  ?? Oct  2015   Amanda Balderrama   amandagbalderrama@gmail.com     3.2
%  12 Nov  2015   Amanda Balderrama   amandagbalderrama@gmail.com     4
%  20 Nov  2015   Amanda Balderrama   amandagbalderrama@gmail.com     4.1
%     Attempting a more principled method for combined map
%  25 Jan  2016   Amanda Balderrama   amandagbalderrama@gmail.com     5
%     Skipped development on v4+ since version didn't seem to be
%     functioning properly
%  25 Feb  2016   Amanda Balderrama   amandagbalderrama@gmail.com     5.1
%     Script debugging for LacSac data (debugs @ 203, 211, 340)
%     -- COMMENTED TEXT CLEANED OUT in v5.2
%  28 Feb  2016   Amanda Balderrama   amandagbalderrama@gmail.com     5.2
%  27 May  2016   Amanda Balderrama   amandagbalderrama@gmail.com     6
%     New attempt to incorporate Kmeans for segmentation
%  12 Jan  2017   Amanda Balderrama   amandagbalderrama@gmail.com     6.1
%     Fixing bugs

if ~exist('hdrtxt')
  hdrtxt = [];
end

[d,hdrtxt,d_div,line_types,type_rows] = divMIMSstruct(d,hdrtxt);
%---> Define data relevant parameters
elnames = fieldnames(rmfield(d,skip_fields));
nel = length(elnames);
delcell = struct2cell(rmfield(d,skip_fields));
linelabels = hdrtxt(2:end,1);
uselines = find(cellfun(@(x) ischar(x),linelabels));
linelabels = linelabels(uselines);
labellen = cellfun(@(x) length(x),linelabels);
lineCO = strfind(lower(linelabels),'_line');
lineCOempt = cellfun(@(x) isempty(x),lineCO);
lineCO(lineCOempt) = num2cell(labellen(lineCOempt) + 1);
relhtxt = cellfun(@(x,y) x(1:(y-1)),linelabels,lineCO,'uniformoutput',0);
nisti = cellfun(@(x) isempty(x),relhtxt);
relhtxt(nisti) = {'NIST612'};
[line_types,ia,ic] = unique(relhtxt,'stable');
ntype = length(ia);
% typerngs is a matrix with (# blocks) rows x 4 columns. A block of data is
% a sequence of lines that all have the same line_type label for a total of
% ntype unique line labels. The columns are defined as follows:
%   {'Block start row','Block end row','# rows in block','Block type label'};
%  'Block type label' is the index in the variable line_types which
%  summarizes the label of all the lines in a given block (for instance:
%  line_types{typerngs(BLOCK_INDEX,4)} will give the label of the data in
%  the block BLOCK_INDEX
typerngs = uselines( find(diff(ic)~=0));
typerngs = [[1;typerngs(1:end)+1],[typerngs;length(hdrtxt(2:end,1))]];
typerngs = [typerngs,typerngs(:,2) - typerngs(:,1) + 1, ic(typerngs(:,2))];
colstart = 1; colend = 2; colLen = 3; colLabel = 4;

%---> Initialize workspace
[M,N] = size(delcell{1});
ncol = 4;
dperblkmask = zeros(M,N);
dtypelabel = zeros(M,N);
dchMask = nan([M,N,ncol]);
ntypeblks = size(typerngs,1);
firstfivesec = find(nanmean(d.Time(:,1:round(N/2))) <= 5);

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

if strmatch('rmout',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  laoffel = PropertyVal{strmatch('rmout',PropertyNames)};
else
  nistel = [];
  otherel = [];
  laoffel = cell(ntype,1);
  nistline = cellfun(@(x) ~isempty(strfind(x,'NIST')),line_types);
  laoffel(nistline) = {nistel};
  laoffel(nistline == 0) = {otherel};
end

if strmatch('tsuel',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  tsueli = PropertyVal{strmatch('tsuel',PropertyNames)};
else
  tsueli = 'S_34';
end

if ~ismember(tsueli,elnames)
  tsueli = 'S_33';
end

if strmatch('label',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  labelmethod = PropertyVal{strmatch('label',PropertyNames)};
else
  labelmethod = 'ch';
end

if strmatch('ncc',PropertyNames)
  ncc = PropertyVal{strmatch('ncc',PropertyNames)};
else
  ncc = [];
end

if strmatch('runv',PropertyNames)
  runv = PropertyVal{strmatch('runv',PropertyNames)};
else
  runv = 0;
end

tissueMask = cellfun(@(x) false(size(x.Time)),d_div,'uniformoutput',0);
laONmask = cellfun(@(x) false(size(x.Time)),d_div,'uniformoutput',0);
labelMask = cellfun(@(x) false(size(x.Time)),d_div,'uniformoutput',0);
ncol = 4;
[M,N] = size(d.Time);
dchMask = nan(M,N,ncol);

if plotimgs
  spcols = ceil( ntype/2 );
  fh = figure('name',strcat(d.dataset,' Laser ON Mask'),...
    'position',[50,50,350*spcols,475*2]);
else
  fh = [];
end

if runv == 2
  [d,hdrtxt,d_div,line_types,type_rows,labelMask,laONmask,tsumask,dchMask] = ...
    labelMIMSstruct_v2(d,hdrtxt,varargin{:});
  for t = 1:ntype
    Nt = size(labelMask{t},2);
    %dperblkmask(type_rows{t},1:Nt) = labelMask
  end
else
  for t = 1:ntype
    %---> Identify relevant parameters in datasetgiven type label
    %disp(sprintf('Type loop: t = %d',t))
    alltyperows = [];
    blkrows = find(typerngs(:,colLabel) == t);
    nblks = length(blkrows);
    typelabel = line_types{typerngs(blkrows(1),colLabel)};
    %---> For each type label, generate a mask
    if ~isempty(strfind(lower(typelabel),'calib'))
      % calib type: for each block of a given calibration segment, create a
      % mask
      for b = 1:nblks
        %disp(sprintf('block loop: b = %d',b))
        drows = typerngs(blkrows(b),colstart):typerngs(blkrows(b),colend);
        alltyperows = [alltyperows,drows];
        [tmask,keepcol,l10muimg,laONeli,fgeli] = makela0FGmask(drows,delcell,firstfivesec,0);
        blkr = typerngs(blkrows(b),1):typerngs(blkrows(b),2);
        blkmask = tmask;
        dperblkmask(blkr,:) = blkmask;
        dtypelabel(blkr,keepcol) = typerngs(blkrows(b),colLabel);
      end
    elseif ~isempty(strfind(lower(typelabel),'nist'))
      drows = [];
      for b = 1:nblks
        drows = [drows,typerngs(blkrows(b),colstart):typerngs(blkrows(b),colend)];
      end
      alltyperows = drows;
      [tmask,keepcol,l10muimg,laONeli,fgeli] = makela0FGmask(drows,delcell,firstfivesec,0);
      for b = 1:nblks
        blkr = typerngs(blkrows(b),1):typerngs(blkrows(b),2);
        blkmask = tmask((sum(typerngs(blkrows(1:b-1),3))+1):sum(typerngs(blkrows(1:b),3)),:);
        dperblkmask(blkr,:) = blkmask;
        dtypelabel(blkr,keepcol) = typerngs(blkrows(b),colLabel);
      end
      %     dblk = cellfun(@(x) x(:,keepcol),normtypedat, 'uniformoutput',0);
      %     dblkv = cellfun(@(x) x(:),dblk,'uniformoutput',0);
      %     dblkm = cell2mat(dblkv');
      %     kmout = kmeans(dblkm,2);
      %     figure; imagesc(reshape(kmout,nr,Np));
    elseif ~isempty(strfind(lower(typelabel),'missingfile'))
      tmask = nan(1);
      drows = [];
      for b = 1:nblks
        drows = [drows,typerngs(blkrows(b),colstart):typerngs(blkrows(b),colend)];
      end
      alltyperows = drows;
      blkr = drows;
      blkmask = tmask;
      dperblkmask(blkr,:) = blkmask;
      dtypelabel(blkr,keepcol) = typerngs(blkrows(b),colLabel);
    else
      drows = [];
      for b = 1:nblks
        drows = [drows,typerngs(blkrows(b),colstart):typerngs(blkrows(b),colend)];
      end
      alltyperows = drows;
      [tmask,keepcol,l10muimg,laONeli,fgeli] = makela0FGmask(drows,delcell,firstfivesec,1);
      nblks = 1;
      blkr = drows;
      blkmask = tmask;
      dperblkmask(blkr,:) = blkmask;
      dtypelabel(blkr,keepcol) = typerngs(blkrows(b),colLabel);
    end
    type_rows{t} = alltyperows(:);
    tissueMask{t} = dperblkmask(alltyperows,keepcol) == 4;
    laONmask{t} = dperblkmask(alltyperows,keepcol) > 2;
    labelMask{t} = dperblkmask(alltyperows,keepcol);
  end %---> end type loop
  
  for c = 1:ncol
    nanmat = nan(M,N);
    nanmat(dperblkmask == c) = 1;
    dchMask(:,:,c) = nanmat;
  end
end

outStruct.line_types = line_types;
outStruct.type_rows = type_rows;
outStruct.elnames = elnames;
outStruct.labelblocks = labelMask;
outStruct.typerowhdrs = {'Block start row','Block end row','# rows in block','Block type label'};
outStruct.typerowInfo = typerngs;
outStruct.colstart = 1;
outStruct.colend = 2;
outStruct.colLen = 3;
outStruct.colLabel = 4;
outStruct.dLabelMask = dperblkmask;
outStruct.dTypeImg = dtypelabel;
outStruct.laONmask = laONmask;
outStruct.tissueMask = tissueMask;
outStruct.dchMask = dchMask;
outStruct.dataset = d.dataset;

varargout{1} = 'ERROR: use v2, previously laONmask';
varargout{2} = 'ERROR: use v2, previously tsumask';
varargout{3} = 'ERROR: use v2, previously dchmask';

%==================================================================
% SUBFUNCTIONS
%==================================================================
  function [tmask,keepcol,l10muimg,laONeli,fgeli] = makela0FGmask(drows,datcell,minla0col,TRImodetype)
    targetsig = 1;
    targetmu = 0;
    numel = length(datcell);
    
    normtypedat = cellfun(@(x) x(drows,:)./nanmax(nanmax(x(drows,:))),datcell,'uniformoutput',0);
    mutypeimg = nanmean(cat(3,normtypedat{:}),3);
    dsize = size(mutypeimg);
    smallestcc = prod(dsize)*0.05^2;
    l10muimg = mutypeimg;
    l10muimg(mutypeimg == 0) = nanmin(mutypeimg(mutypeimg ~=0));
    l10muimg = log10(l10muimg);
    noNANmuimg = ~isnan(l10muimg);
    noNANmuimg = noNANmuimg .* mutypeimg > 0;
    keepcol = find(sum(noNANmuimg) > 0); % identifies the # of columns in the subtype image
    typedat = cellfun(@(x) x(drows,keepcol),datcell,'uniformoutput',0);
    onrowmin = sum(noNANmuimg);
    onrowmin = min(onrowmin(onrowmin > dsize(1)*.2));
    Mp = length(drows); 
    Np = length(keepcol);
    
    [kmout,Ckm] = kmeans(l10muimg(noNANmuimg),2,'emptyaction','singleton','replicates',3);
    vmaskKM = zeros(Mp,Np);
    vmaskKM(noNANmuimg) = kmout;
    
    if TRImodetype
      nmodes = 3;
    else
      nmodes = 2;
    end
    
    [vmask,bimodeOut] = bimodalHistMask(l10muimg,'plotimgs',0,'ncc',1,'nmodes',nmodes,'minpnt',1);
    
    kmout = kmeans(l10muimg(noNANmuimg),2,'emptyaction','singleton','replicates',3);
    vmaskKM = ones(Mp,Np);
    vmaskKM(noNANmuimg) = kmout;
    if median(median(vmaskKM(:,minla0col))) == 2
      vmaskKM(vmaskKM == 2) = 0;
    else
      vmaskKM = vmaskKM - 1;
    end
    vmaskpad = ones(size(vmaskKM,1)+2,size(vmaskKM,2));
    vmaskpad(2:end-1,:) = vmaskKM;
    vmaskpad = imfill(vmaskpad,'holes');
    vmaskKM = vmaskpad(2:end-1,:);
    vmaskKM = [vmaskKM,zeros(Mp,dsize(2)-Np)];
    
    vmaskpad = ones(size(vmask,1) + 2,size(vmask,2));
    vmaskpad(2:(end-1),:) = vmask;
    vmaskpad = imfill(vmaskpad,'holes');
    vmask = vmaskpad(2:(end-1),:);
    svmask = sum(vmask); %<-- sum of all rows: # "on" points per column
    imgarea = Mp*Np;
    fgeli = [];
    laONeli = [];
    vmask2 = vmask;
    strmod = ''; cbh = [];
    %close all
    %figure; subplot(1,3,1); imagesc(l10muimg(:,keepcol)); subplot(1,3,2); imagesc(vmask(:,keepcol)); subplot(1,3,3); plot(bimodeOut.histbins,bimodeOut.histcnts,'marker','o'); hold on; axis tight; plot(bimodeOut.peakI,bimodeOut.peakhistval,'pr'); plot(bimodeOut.lowerCO*[1,1],[0,max(bimodeOut.peakhistval)],':k');
    %figure('position',[24   501   924   496]); subplot(2,3,1); imagesc(l10muimg(:,keepcol));
    %subplot(2,3,2);plot(bimodeOut.histbins,bimodeOut.histcnts,'marker','o'); hold on; axis tight; plot(bimodeOut.peakI,bimodeOut.peakhistval,'pr'); plot(bimodeOut.lowerCO*[1,1],[0,max(bimodeOut.peakhistval)],':k');
    %subplot(2,3,4); imagesc(vmaskKM(1:Mp,1:Np)); title('K-Means mask');
    %subplot(2,3,5); imagesc(vmask(1:Mp,1:Np)); title('bimodal hist mask');
    if sum((svmask > onrowmin*.8)) < Np*0.5 || ...
        sum(svmask(:,minla0col)) > length(minla0col)*onrowmin*0.5 || ... % says at least 50% of the minimum laser off area indicated should be labeled as 0's, otherwise the mask isn't representative 
        TRImodetype == 1
      % Masking technique used didn't yield regions with high "block
      % area" membership and another masking scheme should be used
      nmodes = 2;
      laONall = zeros([Mp,Np,numel]);
      for el = 1:length(datcell)
        %disp(sprintf('makela0FGmask: el = %d',el));
        dmaskraw = datcell{el}(drows,keepcol);
        dmaskraw(isnan(dmaskraw)) = nanmin(dmaskraw(:));
        [~,~,d4mask] = auto_thresh(dmaskraw,'auto',[]);
        if ~any(d4mask(:))
          d4mask = dmaskraw;
        end
        Ad = sqrt(targetsig/nanstd(d4mask(:))^2);
        Bd = targetmu - Ad*nanmean(d4mask(:));
        typedat{el} = d4mask.*Ad + Bd;
        [elm2,elI2] = bimodalHistMask(d4mask,'Nmedfilter',1,'quantfrac',0,'plotimgs',0,...
          'nbins',50,'nmodes',2,'fill',0,'close',0);%,'plotimgs',plotimgs,'ncc',ncc 
        %set(gcf,'position',[452   458   434   333],'name',['ElNum',num2str(el)]);
        %set(gca,'fontsize',16); grid on;
        %figure; subplot(1,3,1); imagesc(d4mask); subplot(1,3,2); imagesc(elm2); subplot(1,3,3); plot(elI2.histbins,elI2.histcnts,'marker','o'); hold on; axis tight; plot(elI2.peakI,elI2.peakhistval,'pr'); plot(elI2.lowerCO*[1,1],[0,max(elI2.peakhistval)],':k');
        if length(elI2.histbins) > 6 && ...
            ~(elI2.peakI(1) == elI2.histbins(1) && elI2.peakI(2) == elI2.histbins(end))
          [elm3,elI3] = bimodalHistMask(d4mask,'Nmedfilter',1,'quantfrac',0,'plotimgs',0,...
            'nbins',50,'nmodes',3,'fill',0,'close',0);%,'plotimgs',plotimgs,'ncc',ncc   
          %figure; subplot(1,3,1); imagesc(d4mask); subplot(1,3,2); imagesc(elm3); subplot(1,3,3); plot(elI3.histbins,elI3.histcnts,'marker','o'); hold on; axis tight; plot(elI3.peakI,elI3.peakhistval,'pr'); plot(elI3.lowerCO*[1,1],[0,max(elI3.peakhistval)],':k');
        else
          elm3 = elm2;
        end
        
        ed = abs(elm2 - elm3);
        
        if sum(sum(imopen(ed,strel('disk',2))))/imgarea > 0.05 && sum(elm3(:)) <= (imgarea - length(minla0col)*Mp)% this indicates that the mask is picking out some "foreground" feature which is occupying at least 25% of the imaging area
          laONall(:,:,el) = elm3;
          elbimode = elI3;
        else
          laONall(:,:,el) = elm2;
          elbimode = elI2;
        end
        
        colsum = sum(laONall(:,:,el));
        elcc = bwconncomp(laONall(:,:,el));
        complen = cellfun(@(x) length(x),elcc.PixelIdxList);
        if isempty(complen)
          complen = 0;
        end
        %[~,cent,sdel] = kmeans(d4mask(:),elm3(:));
        %m1 = nanmean(d4mask(elm3 == 1));
        %s1 = nanstd(d4mask(elm3 == 1));
        %m0 = nanmean(d4mask(elm3 == 0));
        %s0 = nanstd(d4mask(elm3 == 0));
        %figure; subplot(1,3,1); imagesc(d4mask(:,keepcol)); subplot(1,3,2); imagesc(laONall(:,:,el)); subplot(1,3,3); plot(elbimode.histbins,elbimode.histcnts,'marker','o'); hold on; axis tight; plot(elbimode.peakI,elbimode.peakhistval,'pr'); plot(elI2.lowerCO*[1,1],[0,max(elbimode.peakhistval)],':k');
        if sum((colsum > Mp*.8))>Np*0.5 && sum(colsum(minla0col)) < length(minla0col)*onrowmin*0.5
          laONeli = [laONeli,el];
          %imoverlay(laONall(:,:,el),'',laONall(:,:,el),datcell{el}(drows,keepcol),'','npanels',1,'displaymode','rgbgray','colormap','gray','transparency',0.3)
          %figure('position',[817   718   423   260],'name',['ElNum',num2str(el),'_Bkgnd']); imagesc(laONall(:,:,el)); colormap(gray); 
          %figure('position',[817   718   423   260],'name',['ElNum',num2str(el),'_BkgndImg']); imagesc(datcell{el}(drows,keepcol)); hold on; contour(laONall(:,:,el),'r'); 
          %set(gcf,'position',[817   718   423   260],'name',['ElNum',num2str(el),'_Bkgnd']);set(gca,'position',[0,0,1,1]);
        elseif max(complen) >= imgarea*0.25^2 && sum(sum(laONall(:,minla0col,el))) <= length(minla0col)*Mp*.5%sum(sum(laONall(:,:,el))) <= (imgarea - length(minla0col)*Mp)% this indicates that the mask is picking out some "foreground" feature which is occupying at least 25% of the imaging area
          %figure; subplot(1,3,1); imagesc(d4mask(:,keepcol)); subplot(1,3,2); imagesc(laONall(:,:,el)); subplot(1,3,3); plot(elbimode.histbins,elbimode.histcnts,'marker','o'); hold on; axis tight; plot(elbimode.peakI,elbimode.peakhistval,'pr'); plot(elI2.lowerCO*[1,1],[0,max(elbimode.peakhistval)],':k');
          fgeli = [fgeli,el];
          %figure('position',[817   718   423   260],'name',['ElNum',num2str(el),'_Tissue']); imagesc(laONall(:,:,el)); colormap(gray); set(gca,'position',[0,0,1,1]);
          %figure('position',[817   718   423   260],'name',['ElNum',num2str(el),'_TissueImg']); imagesc(datcell{el}(drows,keepcol)); hold on; contour(laONall(:,:,el),'r'); 
          %set(gca,'position',[0,0,1,1]); set(gcf,'position',[817   718   423   260],'name',['ElNum',num2str(el),'_Tissue']);
        end
      end
      
      if ~isempty(laONeli)
        lsum = sum(laONall(:,:,laONeli),3); %assumes the first element is carbon which is a negative of the mask
        coval = round(mean(mean(lsum(:,minla0col))));% mode(lsum(lsum>min([1,length(laONeli)-1]))) - 1;
      else
        lsum = sum(laONall,3); %assumes the first element is carbon which is a negative of the mask
        coval = round(mean(mean(lsum(:,minla0col))));%mode(lsum(lsum>min(lsum(:)))) - 1;
      end
      coval = quantile(reshape(lsum(:,minla0col),[prod(size(lsum(:,minla0col))),1]),0.75);
      
      vmask2 = reshape(kmeans(lsum(:),2,'emptyaction','singleton','replicates',3),size(lsum));
      if median(median(vmask2(:,minla0col))) == 2
       vmask2(vmask2 == 2) = 0;
      else
       vmask2 = vmask2 - 1;
      end
 
      vmask2pad = ones(size(vmask2,1)+2,size(vmask2,2));
      vmask2pad(2:end-1,:) = vmask2; 
      vmask2pad = imfill(vmask2pad,'holes');
      vmask2 = vmask2pad(2:end-1,:);
      vmask2 = [vmask2,zeros(size(vmask2,1),size(vmask,2)-size(vmask2,2))];

      %figure; subplot(1,2,1); imagesc(lsum(:,keepcol)); subplot(1,2,2); imagesc(vmask2(:,keepcol));
      if abs(sum(vmask(:)-vmask2(:))) > imgarea*0.05 &&... 
          sum(vmask2(:)) > sum(vmask(:)) 
        % difference between the area indicated by the two masks is greater
        % than 5% of the image area AND the mask 2 method has larger area
       vmask = vmask2;
      elseif sum(sum(vmask2) >= onrowmin) > sum(svmask >= onrowmin) % This covers the case where the original vmask has gaps 
        % there are more laser on rows in mask 2 than in mask 1 (ie mask 2
        % is more rectangular and blocky)
        vmask = vmask2;
      elseif abs(sum(vmask(:)-vmask2(:))) < imgarea*0.1 &&...
          abs(sum(sum(vmask2(:,minla0col))) - sum(svmask(minla0col))) > length(minla0col)*dsize(1)*0.05
        vmask = vmask2;
      elseif sum(svmask(:,minla0col)) > length(minla0col)*onrowmin*0.5
        vmask = vmask2;
      end
      %vmask = vmask2;
      svmask = sum(vmask);
      %laoffval = mode(mode(lsum(:,1:fiveseccol)));
    end
    
    svmask = svmask(:);
    objcol = kmeans(svmask(1:Np),2,'emptyaction','singleton','replicates',3) - 1;
    objcol = find(objcol ~= objcol(1));
    % --> old, replaced 26 feb 2016: objcol = find(svmask >= onrowmin*.8);% max(svmask));
    firsthalfi = round(mean(objcol));
    
    tmask = zeros(size(mutypeimg));
    tmask(:,1:firsthalfi) = abs(vmask(:,1:firsthalfi) - 1) + tmask(:,1:firsthalfi);
    tmask = tmask + vmask.*3;
    tmask(tmask == 0) = 2;
    tmask(~noNANmuimg) = 0;
    [la0_1r,la0_1c] = find(tmask == 1);
    [la0_2r,la0_2c] = find(tmask == 2);
    [bgr,bgc] = find(tmask == 3);
    
    %if quantile(la0_1c,0.75) > quantile(bgc,0.05)
    %  tmask(:,1:firsthalfi) = 3;
    %  tmask(:,1:max([(min(bgc)-1),1])) = 1;
    %end
    
    %if quantile(la0_2c,0.25) < quantile(bgc,0.95)
    %  tmask(:,(firsthalfi + 1):(Np-length(minla0col))) = 3;
    %  tmask(:,(Np - length(minla0col)):Np) = 2;
    %end
    
    %<< FOR DEBUGGING, UNCOMMENT TO PLOT OUTPUTS (also uncomment line 388)>>
    % figure; subplot(3,2,1); imagesc(l10muimg(:,keepcol)); colorbar; subplot(3,2,2); hist(l10muimg(~isnan(l10muimg)),100); hold on; plot(bimodeOut.lowerCO,0,'rp'); subplot(3,2,3); imagesc(vmask(:,keepcol)); subplot(3,2,4); imagesc(vmask2(:,keepcol)); subplot(3,2,5); imagesc(tmask(:,keepcol)); 
    
    if TRImodetype
      imgarea = sum(tmask(:) == 3);  % added 20160404 
      modmask = zeros(dsize);
      la0med = nanmedian(l10muimg(tmask <= 2));
      la0std = nanstd(l10muimg(tmask <= 2));
      if abs(bimodeOut.peakI(1) - la0med) > la0std
        % indicates that the first peaks detected in tri-modal detection
        % doesn't correspond to the laser off period
        [modmask,bimodeOut2] = bimodalHistMask(l10muimg,'nmodes',2);
      else
        % indicates that the first peak detected corresponds to laser off and
        % the others likely correspond to background and tissue
        modmask(l10muimg > mean(bimodeOut.peakI(2:end))) = 1;
      end
      ccfg = bwconncomp(modmask);
      lenccfg = cellfun(@(x) length(x),ccfg.PixelIdxList);
      keepcc = find(lenccfg >= smallestcc);
      modmask = zeros(size(mutypeimg));
      modmask(cat(1,ccfg.PixelIdxList{keepcc})) = 1;
      %figure; imagesc(imclose(fgmask,strel('disk',1)));%[max([round(dsize(1)*.01),3]),min([round(dsize(2)*.01),5])])));
     
      if isempty(fgeli)
        lsum = sum(laONall,3);
      else
        lsum = sum(laONall(:,:,fgeli),3);
      end
      %figure; bar(elbimode.histbins, elbimode.histcnts)
      %figure; subplot(2,2,1); imagesc(l10muimg(:,keepcol)); subplot(2,2,2); imagesc(modmask(:,keepcol)); subplot(2,2,3); imagesc(lsum); subplot(2,2,4); imagesc(fgmask(:,keepcol));
            
      [fgmask,elbimode] = bimodalHistMask(lsum,'Nmedfilter',1,'quantfrac',0,'plotimgs',0,...
        'nbins',length(fgeli)+1,'nmodes',2,'fill',0,'close',0);
      %tsumat = cell2mat(cellfun(@(x) x(:),typedat(fgeli),'uniformoutput',0)');
      %ktsu = kmeans(tsumat,2);
      %figure('position',[680   565   838   413]); subplot(2,3,1); imagesc(modmask(:,keepcol)); title('bimodalHist(mean log image)'); 
      %subplot(2,3,2); imagesc(reshape(ktsu,Mp,Np)); title('kmeans(FG elements)');
      %subplot(2,3,3); imagesc(fgmask(:,keepcol)); title('bimodalHist(FG mask)');
      %subplot(2,3,4); imagesc(l10muimg(:,keepcol)); subplot(2,3,5); imagesc(reshape(mean(tsumat,2),Mp,Np)); subplot(2,3,6); imagesc(lsum(:,keepcol));
      % --> legacy method for determining cutoff threshold
      %coval = mode(lsum(lsum>min([1,length(fgeli)-1]))) - 1;%quantile(lsum(tmask <= 2),0.9) + 1;%
      %fgmask(:,keepcol) = lsum >= coval;
      %bgmod = quantile(lsum(fgmask == 0),0.9);
      %fgmask(:,keepcol) = lsum >= (coval+bgmod)/2;
      fgmask = lsum >= round(elbimode.lowerCO);
      dFM = (2*sum(sum(fgmask(:,keepcol).*modmask(:,keepcol))))/(sum(modmask(1:Np*Mp)) + sum(fgmask(1:Np*Mp)));
      dTM = (2*sum(sum(modmask(:,keepcol).*(tmask(:,keepcol)==3))))/(sum(sum(tmask(:,keepcol) == 3)) + sum(modmask(1:Np*Mp)));
      dFT = (2*sum(sum(fgmask(:,keepcol).*(tmask(:,keepcol)==3))))/(sum(sum(tmask(:,keepcol) == 3)) + sum(fgmask(1:Np*Mp)));
      
      % Requires tuning: Objective is to distinguish use of modmask in the
      % case 
      if abs(sum(fgmask(:)) - sum(modmask(:))) > imgarea*0.1 && ...
          sum(modmask(:)) > imgarea*0.1 && sum(modmask(:)) < imgarea*0.9 
        if dFT < 0.9 
          fgmask = modmask;
        end
      elseif ~any(fgmask(:))
        fgmask = modmask;
      end
      ccfg = bwconncomp(fgmask);
      lenccfg = cellfun(@(x) length(x),ccfg.PixelIdxList);
      keepcc = find(lenccfg >= smallestcc);
      fgmask = zeros(size(mutypeimg));
      fgmask(cat(1,ccfg.PixelIdxList{keepcc})) = 1;
      
      cc0 = bwconncomp(abs(fgmask - 1));
      lencc0 = cellfun(@(x) length(x),cc0.PixelIdxList);
      fillcc = find(lencc0 <= smallestcc*0.1);
      fgmask(cat(1,cc0.PixelIdxList{fillcc})) = 1;
      [fgr,fgc] = find(fgmask == 1);
      
      %[h_l10,p_l10,ci_l10,stats_l10] = ttest2(l10muimg(modmask == 1),l10muimg(modmask(:,keepcol) == 0),'Vartype','unequal');
      %[h_lsum,p_lsum,ci_lsum,stats_lsum] = ttest2(l10muimg(fgmask == 1),l10muimg(fgmask(:,keepcol) == 0),'Vartype','unequal');
      
      
      %fgmask = imclose(fgmask,strel('disk',1));
      %figure; imagesc(imclose(fgmask,strel('disk',1)));%[max([round(dsize(1)*.01),3]),min([round(dsize(2)*.01),5])])));
      if quantile(la0_1c,0.75) > quantile(bgc,0.05)
        tmask(:,1:firsthalfi) = 3;
        tmask(:,1:max([(min(fgc)-1),1])) = 1;
      end
      
      if quantile(la0_2c,0.25) < quantile(bgc,0.95) ||...
        sum(la0_2c < quantile(bgc,0.99))/length(la0_2c) > 0.1 % more than 10% of the columns are "laser off" BEFORE the end of the tissue, indicates and erro
        tmask(:,(firsthalfi + 1):Np) = 3;
        tmask(:,max([(Np - length(minla0col)),max(fgc)]):Np) = 2;
      end
      
      tmask(fgmask & tmask == 3) = 4;
      %<< FOR DEBUGGING, UNCOMMENT TO PLOT OUTPUTS>> \
      %subplot(3,2,6); imagesc(tmask(:,keepcol))
    end
  end

end
