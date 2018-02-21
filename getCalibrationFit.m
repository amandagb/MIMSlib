function [Calstruct,dOut,varargout] = getCalibrationFit(d_div,line_types,laONmask,type_rows,varargin)
%% function: getCalibrationFit
% Description: **IMPORTANT NOTE** No outlier removal is done to compute
% means and medians of data. Indeed, the median should in theory be the
% statistic which is more robust to outliers in the data.
% 
% Example:
%     >> [cals,dcal] = getCalibrationFit(d_div,line_types,laONmask);
%     >> [rawcal,drawcal] = getCalibrationFit(d_div,line_types,labelInfo.laONmask,type_rows,'plotimgs',1,'calibstr','Delta_calib');
%
%
% INPUTS ----------------------------------------------------------------
% datastruct:    MIMS structure OR a string indiciating the dataset the
% user would like to use.
% linelabels:
% varargin - 'PropertyName','PropertyValue'
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'plotbounds': logical indicating whether to plot bounds or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'ElLimits': nel x 2 matrix with rows corresponding to the
%   [min(el),max(el)] values for which the ylimits of the calibration
%   curves should be [DEFAULT = []]
%   • 'sqmask': logical indicating whether the calibration data should be
%   masked with a square [DEFAULT = 0]
%   • 'calcol': vector indicating the columns to include in the calibration
%   analysis [DEFAULT = trues in Mask]
%   • 'normbyel': string indicating isotope to normalize by [DEFAULT = []]
%   • 'calibstr': string indicating the "header string" to search for
%   [DEFAULT = 'Calib']
%
% OUTPUTS ---------------------------------------------------------------
% Calstruct: Stucture with relevant output data
%   • Elnames:
%   • LineSlopeInt_Means: nel x 2 matrix with rows corresponding to the
%   [slope,intercept] values for the element corresponding to that row
%   • LineFitR2_Means = rsqmean;
%   • LineSlopeInt = p;
%   • LineFitR2 = rsq;
% dOut:
%   • d_div: L x 1 cell of structures where each structure is the
%        basic MIMS structure associated with the label
%   • line_types
%   • type_rows
%   • CalMasks
%   • dinMasks
%   • ppmVal
%   • CalLabels
%   • dmeans
%
%  Date           Author              E-mail                      Version
%  24 Mar  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
%   1 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%     Separated portions of the script into multiple functions since each
%     individually is commonly used
%  30 Sept 2015   Amanda Balderrama   amandagbalderrama@gmail.com     3
%     Incorporated processing for NIST sections as well
%  1  Oct  2015   Amanda Balderrama   amandagbalderrama@gmail.com     3.1
%     Additional varargin which allows for normalization by a specific
%     isotope -- experimental for protein standards
% 18  Feb  2016   Amanda Balderrama   amandagbalderrama@gmail.com     4
%     Added the ability to specify the calibration label to look for in the
%     data

compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 1;
end

if strmatch('plotbounds',PropertyNames)
  bnds = PropertyVal{strmatch('plotbounds',PropertyNames)};
else
  bnds = 0;
end

if strmatch('plotmean',PropertyNames)
  plotmean = PropertyVal{strmatch('plotmean',PropertyNames)};
else
  plotmean = 1;
end

if strmatch('plotmed',PropertyNames)
  plotmed = PropertyVal{strmatch('plotmed',PropertyNames)};
else
  plotmed = 1;
end

if strmatch('plotpnts',PropertyNames)
  plotpnts = PropertyVal{strmatch('plotpnts',PropertyNames)};
else
  plotpnts = 1;
end

if strmatch('saveimgs',PropertyNames)
  saveimgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  saveimgs = 0;
end

if strmatch('plotel',PropertyNames)
  plotel = interp_data_elemrng(d_div{1},PropertyVal{strmatch('plotel',PropertyNames)});
else
  plotel = interp_data_elemrng(d_div{1},[]);
end

if strmatch('xlabelstr',PropertyNames)
  xlab = PropertyVal{strmatch('xlabelstr',PropertyNames)};
else
  xlab = 'C = Standard Concentration [ppm]';
end

if strmatch('ylabelstr',PropertyNames)
  ylab = PropertyVal{strmatch('xlabelstr',PropertyNames)};
else
  ylab = 'I = Mean Intensity [cps]';
end

if strmatch('ellim',PropertyNames)
  ellims = PropertyVal{strmatch('ellim',PropertyNames)};
else
  ellims = [];
  %   ellims = [4e5,2e6;...'C_13'
  %     0e5,10e5;...'Ce140'
  %     0,1e6;...'Cu63'
  %     0,2e5;...'Cu65'
  %     0,6e5;...'Eu151'
  %     0,1e4;...'Fe57'
  %     0,1.5e5;...'Gd157'
  %     2.5e5,7.5e5;...'K_39'
  %     0,7e5;...'La138'
  %     0,3e5;...'Mn55'
  %     0,6e4;...'Mo95'
  %     1e5,7e5;...'P_31'
  %     0.5e4,2.5e4;...'S_33'
  %     2e4,9e4;...'S_34'
  %     3.5e7,5.5e7;...'S_36'
  %     0,5e5;...'Sr88'
  %     0,2e5;...'Ti48'
  %     0,6e5;...'Y_89'
  %     2e4,7e4;...'Zn64'
  %     0,2.25e4;...'Zn66'
  %   ];
end

if strmatch('sqmask',PropertyNames)
  sqmask = PropertyVal{strmatch('sqmask',PropertyNames)};
else
  sqmask = 0;
end

if strmatch('calcol',PropertyNames)
  calcol = PropertyVal{strmatch('calcol',PropertyNames)};
else
  calcol = [];
end

if strmatch('normbyel',PropertyNames)
  normbyel = PropertyVal{strmatch('normbyel',PropertyNames)};
else
  normbyel = [];
end

if strmatch('calibstr',PropertyNames)
  calibstr = PropertyVal{strmatch('calibstr',PropertyNames)};
else
  calibstr = [];
end
%% Reads raw data
if isstr(d_div)
  [~,~,d_div,line_types,~,~,laONmask] = labelMIMSstruct(d_div);
end

delonly = rmfield(d_div{1},skip_fields);
elnames = sort(fieldnames( delonly ));
nel = length(elnames);

%% Uses the calibration mask to extract data for all elements in each channel
calibpos = strfind(lower(line_types),lower(calibstr));
calibcells = find(cellfun(@(x) ~isempty(x),calibpos));
nstnds = length(calibcells);

if sqmask 
  % >>> Create square, equal area masks in all channels
  %   The code in this block looks at the existing foregorund masks
  %   (provided by laONmask variable) and picks the columns to keep as
  %   those that have less than 10% of the lines indicating that the column
  %   is background in that row
  calsize = cell2mat(cellfun(@(x) size(x),laONmask(calibcells),'uniformoutput',0));
  if length(unique(calsize(:,2))) > 1
    %disp('Develop case -- need to resize masks to they''re all the same size')
    maskblk = false([sum(calsize(:,1)),max(calsize(:,2))]);
    startfrom = 0;
    for c = 1:nstnds
      maskblk((sum(calsize(1:c-1,1))+1):sum(calsize(1:c,1)),1:calsize(c,2)) = laONmask{calibcells(c)};
    end
  else
    maskblk = cell2mat(laONmask(calibcells));
  end
  
  if isempty(calcol)
    smb = sum(maskblk); %sum of rows in mask block
    nmr = max(smb); %number of mask block rows
    %smbfg = smb > nmr*0.5; %indices for which more than half the rows indicate that the region is foreground
    fgcol = find(smb >= median(smb));
  else
    fgcol = calcol;
  end
  sqmask = false(size(maskblk)); % <== may require image closing to ensure a single foreground conncented component without holes
  sqmask(:,fgcol) = 1;

  if length(unique(calsize(:,2))) > 1
    for c = 1:nstnds
      laONmask{calibcells(c)} = sqmask((sum(calsize(1:c-1,1))+1):sum(calsize(1:c,1)),1:calsize(c,2));
    end
  else
  laONmask(calibcells) = mat2cell(sqmask,calsize(:,1)',size(maskblk,2));
  end
end

ppmnum = nan(nstnds,1);
caldata = cell(1,nstnds);
for c = 1:nstnds
  calmask = calibcells(c);
  dcal = d_div{calibcells(c)};
  dcalel = struct2cell(rmfield(dcal,skip_fields));
  calnumstr = strsplit(dcal.dataset,'_');
  ppmnum(c) = str2num(calnumstr{end});
  
  npnts = sum(laONmask{calmask}(:));
  if npnts
    caldata{c} = reshape(cell2mat(cellfun(@(x) x(laONmask{calmask}),dcalel,'uniformoutput',0)),...
      npnts,nel);
    
    %if plotimgs --- box plots
    %  figure('position',[284, 530, 1070, 427]);
    %  boxplot(dsmask{c}./repmat(max(dsmask{c}),npnts,1),elnames); hold on;
    %  hold on; plot(mean(dsmask{c})./max(dsmask{c}),'gx');
    %  title(sprintf('%s: Normalized box-whisker & means of data in Standard %d',strrep(dcal.dataset,'_','\_'),ppmnum(c)));
    %end
  else
    caldata{c} = nan(3,nel);
    disp(sprintf('%s does not have any included data in the mask',line_types{calmask}))
  end
end
[ppmnum,reorderrows] = sort(ppmnum);
caldata = caldata(reorderrows);
calInfo.rawdinMasks = caldata;

if ~isempty(normbyel)
  normi = interp_data_elemrng(dcal,normbyel);
  normbyel = elnames{normi};
  caldata = cellfun(@(x) x./repmat(x(:,normi),1,nel),caldata,'uniformoutput',0);
end

dmeans = cell2mat(...
  cellfun(@(x) nanmean(x)',caldata,'uniformoutput',0))';

dmedians = cell2mat(...
  cellfun(@(x) nanmedian(x)',caldata,'uniformoutput',0))';

dq1 = cell2mat(...
  cellfun(@(x) quantile(x,0.25)',caldata,'uniformoutput',0))';
dq3 = cell2mat(...
  cellfun(@(x) quantile(x,0.75)',caldata,'uniformoutput',0))';
diqr = cell2mat(...
  cellfun(@(x) iqr(x)',caldata,'uniformoutput',0))';

dstd = cell2mat(cellfun(@(x) nanstd(x)',caldata,'uniformoutput',0))';

calInfo.dinMasks = caldata;
calInfo.normbyel = normbyel;
calInfo.ppmVal = ppmnum;
calInfo.CalLabels = line_types(calibcells(reorderrows));
calInfo.dmeans = dmeans;
calInfo.dstd = dstd;
calInfo.dmedians = dmedians;
calInfo.dq1 = dq1;
calInfo.dq3 = dq3;
calInfo.diqr = diqr;

%% Linear curve fitting of calibration data
dsmaskmat = cell2mat(caldata(:));
nonNAN = ~isnan(dsmaskmat);
nonNANmeans = ~isnan(dmeans);
x = cell(nstnds,1);
for c = 1:nstnds
  x{c} = repmat(ppmnum(c),size(caldata{c},1),1);
end
x = cell2mat(x);

% p = nan(nel,2);
% rsq = nan(nel,1);
% %S = cell(1,nel);
% %mu = cell(1,nel);
% for i = 1:nel
%   [fobj,gof,outstruct] = fit(x(nonNAN(:,i)),dsmaskmat(nonNAN(:,i),i),'poly1');%'Lower', [-Inf dmeans(1,i)], 'Upper', [Inf dmeans(1,i)]);
%   p(i,:) = [fobj.p1,fobj.p2];
%   rsq(i) = gof.rsquare;
%   %[p(i,:),S{i},mu{i}] = polyfit(x,dsmaskmat(:,i),1);
% end

pmed = nan(nel,2);
rsqmed = nan(nel,1);
%S = cell(1,nel);
%mu = cell(1,nel);
for i = 1:nel
  [fobj,gof,outstruct] = fit(ppmnum(nonNANmeans(:,i)),dmedians(nonNANmeans(:,i),i),'poly1');%'Lower', [-Inf dmeans(1,i)], 'Upper', [Inf dmeans(1,i)]);
  pmed(i,:) = [fobj.p1,fobj.p2];
  rsqmed(i) = gof.rsquare;
  %[p(i,:),S{i},mu{i}] = polyfit(x,dsmaskmat(:,i),1);
end

pmeans = nan(nel,2);
rsqmean = nan(nel,1);
for i = 1:nel
  [fobj,gof,outstruct] = fit(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),'poly1');%'Lower', [-Inf dmeans(1,i)], 'Upper', [Inf dmeans(1,i)]);
  pmeans(i,:) = [fobj.p1,fobj.p2];
  rsqmean(i) = gof.rsquare;
end
uscorei = strfind(d_div{calibcells(1)}.dataset,'_');
datastr = d_div{calibcells(1)}.dataset(1:(uscorei(2)-1));

%% Aggregates data in NIST labeled regions
nisttype = find(cellfun(@(x) ~isempty(x),strfind(line_types,'NIST')));
if ~isempty(nisttype)
  nistrows = type_rows{nisttype};
  nistsections = [1;find(diff(nistrows) > 1)+1];
  nistsections = [nistsections,[nistsections(2:end)-1;length(nistrows)]];
  dnistel = struct2cell(rmfield(d_div{nisttype},skip_fields));
  nnist = sum(laONmask{nisttype}(:));
  maskpntsperrow = [0;sum(laONmask{nisttype},2)];
  nistdata = zeros(sum(maskpntsperrow),nel);
  nnistsec = size(nistsections,1);
  munistsec = zeros(nnistsec,nel);
  mednistsec = zeros(nnistsec,nel);
  signistsec = zeros(nnistsec,nel);
  nistcell = cell(1,nnistsec);
  nistrowsind = zeros(nnistsec,2);
  
  for nr = 1:nnistsec
    maskblk = laONmask{nisttype}(nistsections(nr,1):nistsections(nr,2),:);
    startrow = sum(maskpntsperrow(1:nistsections(nr,1)));
    ndrows = [1:sum(maskpntsperrow((nistsections(nr,1)+1):(nistsections(nr,2)+1)))] + startrow;
    for el = 1:nel
      dnel = dnistel{el}(nistsections(nr,1):nistsections(nr,2),:);
      nistdata(ndrows,el) = dnel(maskblk);
      munistsec(nr,el) = mean(dnel(maskblk));
      mednistsec(nr,el) = median(dnel(maskblk));
      signistsec(nr,el) = std(dnel(maskblk));
    end
    nistcell{nr} = nistdata(ndrows,:);
    nistrowsind(nr,:) = [nistrows(nistsections(nr,1)),nistrows(nistsections(nr,2))];
  end
  
  nistInfo.dinMasks = nistcell;
  nistInfo.dnistrows = nistrowsind;
  nistInfo.dmedians = mednistsec;
  nistInfo.dmeans = munistsec;
  nistInfo.dstd = signistsec;
else
  nistInfo = [];
end

%% Depicting linear fit with raw data points
if plotimgs
  pcol = 4;
  prow = nstnds+length(nisttype);
  spmat = repmat([0:prow-1]',1,pcol).*pcol +repmat(1:pcol,prow,1);
  for i = plotel(:)'%1:nel
    allpnts = [];
    ccel = figure('name',elnames{i},'position',[179    49   986   943]);
    %--------------------------------------------------------
    % Plot for best fit curves
    %--------------------------------------------------------
    set(gcf,'position',[179    49   986   943]);
    spos = spmat(2:prow,1:2);
    subplot(prow,pcol,spos(:)');
    if plotpnts
      plot(x,dsmaskmat(:,i),'.c');
    end
    hold on;
    if bnds
      boundedline(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),dstd(nonNANmeans(:,i),i),'b', 'alpha');
      boundedline(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),dstd(nonNANmeans(:,i),i),'xr', 'alpha');
    else
      if plotmed
        plot(ppmnum,dmedians(:,i),'db');
        medline = pmed(i,1).*ppmnum + pmed(i,2);
        plot(ppmnum,medline,'b');
        allpnts = [allpnts(:);dmedians(:,i);medline];
      end
      if plotmean
        plot(ppmnum,dmeans(:,i),'xr');
        meanline = pmeans(i,1).*ppmnum + pmeans(i,2);
        plot(ppmnum,meanline,'r');
        allpnts = [allpnts(:);dmeans(:,i);meanline];
      end
    end
    
    title({[sprintf('%s',strrep(datastr,'_','\_'))],...
      [sprintf('%s Calibration curves',strrep(elnames{i},'_','\_'))]});
    if isempty(ellims)
      %mxyval = quantile(dsmaskmat(:,i),0.98);
      %ylim([min(dsmaskmat(:,i)),mxyval]);
      ylim([min(allpnts)-eps, max(allpnts)+eps])
    else
      mxyval = ellims(i,2);
      ylim(ellims(i,:));
    end
    cxlim = xlim;
    iylim = ylim;
    %mxyval = ylim(2);

    if ~isempty(nisttype)
      plot(repmat(cxlim,nnistsec,1)',repmat(munistsec(:,i),1,2)','--')
      lcol = lines(nnistsec);
      for nr = 1:nnistsec
        plot((munistsec(nr,i) - pmed(i,2))./pmed(i,1),munistsec(nr,i),'o',...
          'markeredgecolor',lcol(nr,:))
      end
    end
    xlim(cxlim);
    ylim(iylim);
    xlabel(xlab)
    ylabel(ylab)
    cbrng = ylim;
    cbmin = cbrng(1);%max([0,dmeans(1,i)-dstd(1,i)]);
    cbmax = cbrng(2);%min([dmeans(end,i) + 0.5*dstd(end,i),quantile(caldata{end}(:,i),0.75)]);
    cbins = cbmin:(cbmax-cbmin)/24:cbmax;
    
    subplot(prow,pcol,1);
    if plotmed
      text(0,1,...
        sprintf('$$I_{med} = %1.2e C + %1.2e, R^2 = %1.2f$$',pmed(i,1),pmed(i,2),rsqmed(i)),...
        'fontsize',14,'color','b','VerticalAlignment','top','interpreter','latex')
    end
    if plotmean
      text(0,0.5,...
        sprintf('$$I_{mean} = %1.2e C + %1.2e, R^2 = %1.2f$$',pmeans(i,1),pmeans(i,2),rsqmean(i)),...
        'fontsize',14,'color','r','VerticalAlignment','top','interpreter','latex')
    end
    set(gca,'xticklabel','','yticklabel','','color','none','visible','off')
    
    for c = 1:nstnds
      dcal = d_div{calibcells(c)};
      dcalel = struct2cell(rmfield(dcal,skip_fields));
      subplot(prow,pcol,spmat(c,3));
      if ~isempty(normbyel) && i ~= normi
        dcalel{i} = dcalel{i}./dcalel{normi};
      elseif ~isempty(normbyel) && i == normi
        normdraw = cell2mat(calInfo.rawdinMasks(:));
        cbins = auto_thresh(normdraw(:,normi));
        cbmin = cbins(1);
        cbmax = cbins(2);
        cbins = cbmin:(cbmax-cbmin)/24:cbmax;
      end
      
      imagesc(dcalel{i},[cbins(1),cbins(end)]); colormap(gray);%colormap(custom_colormaps('parula'));%
      calmask = laONmask{calibcells(c)};
      [xc,yc] = meshgrid(1:size(calmask,2),1:size(calmask,1));
      xc = xc(calmask == 0);
      yc = yc(calmask == 0);
      npnts = sum(calmask(:));
          
      if c == nstnds
        cbh = colorbar('location','west');
        set(cbh,'YColor','c',...
          'fontsize',12,'fontweight','bold');
      end
      hold on; 
      plot(xc,yc,'.r');%contour(calmask,[0,0],'r');
      title(sprintf('%1.1f ppm',ppmnum(c)));
      
      subplot(prow,pcol,spmat(c,4));
      hc = hist(dcalel{i}(calmask == 1),cbins);
      bar(cbins,hc./npnts);
      xlim([cbmin,cbmax])
      title(sprintf('%1.1f ppm hist',ppmnum(c)));
      xlabel('I [cps]')
    end
    
    if ~isempty(nisttype)
      nistbins = auto_thresh(nistdata(:,i));
      subplot(prow,pcol,spmat(end,3));
      imagesc(dnistel{i},nistbins); hold on;
      [xc,yc] = meshgrid(1:size(laONmask{nisttype},2),1:size(laONmask{nisttype},1));
      xc = xc(laONmask{nisttype} == 0);
      yc = yc(laONmask{nisttype} == 0);
      plot(xc,yc,'.r');%contour(calmask,[0,0],'r');
      title('NIST glass');
      
      nstb = nistbins(1):(nistbins(2)-nistbins(1))/19:nistbins(2);
      hnall = hist(nistdata(:,i),nstb);
      subplot(prow,pcol,spmat(end,4));
      plot(nstb,hnall./sum(hnall),'o');hold on;
      hsec = zeros(nnistsec,length(nstb));
      for nr = 1:nnistsec
        startrow = sum(maskpntsperrow(1:nr));
        ndrows = [1:maskpntsperrow(nr+1)] + startrow;
        hsec(nr,:) = hist(nistdata(ndrows,i),nstb);
        subplot(prow,pcol,spmat(end,3));
        rectangle('position',[0.5,nistsections(nr,1)-.2,size(laONmask{nisttype},2),nistsections(nr,2)-nistsections(nr,1) + .6],...
          'EdgeColor',lcol(nr,:));
      end
      subplot(prow,pcol,spmat(end,4));
      plot(nstb,hsec./repmat(sum(hsec,2),1,length(nstb)));
      xlim(nistbins)
    end
  end
  if saveimgs
    save_open_figures(strcat(compdir.lafldr,datastr,'\Images'),[],[],strcat(datastr,'_CalibrationCurve'))
    %export_fig(strcat('Images\CalibrationCurve_',elnames{i}),...
    %  '-png');
    %close(ccel)
  end
end

% dockf on all
Calstruct.elnames = elnames;
Calstruct.LineSlopeInt_Means = pmeans;
Calstruct.LineFitR2_Means = rsqmean;
Calstruct.LineSlopeInt_Med = pmed;
Calstruct.LineFitR2_Med = rsqmed;
dstr = strrep(d_div{1}.dataset,line_types{1},'');
dOut.dataset = dstr(1:end-1);
dOut.calInfo = calInfo;
dOut.nistInfo = nistInfo;

end
