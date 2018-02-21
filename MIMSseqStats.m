function [dstats,varargout] = MIMSseqStats(d,labelInfo,varargin)
%% function: MIMSseqStats
% Description:
% Example:
%     [dstats,varargout] = MIMSseqStats(d,labelInfo,varargin)
% 
% Required Function: must first run labelMIMSstruct for labelInfo variable
%
% INPUTS ----------------------------------------------------------------
% d:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% hdrtxt: 
% varargin - 'PropertyName','PropertyValue'
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'savestats': logical indicating whether to save stat tables or not [DEFAULT = 0]
%   • 'normbyel':   string indicating the element to normalize by [DEFUALT
%   = none]
%   • 'rmout': logical for removing outliers (1) or not (0) [DEFAULT =  1]
%
% OUTPUTS ---------------------------------------------------------------
% dstats: 
%   • 'ntypeblks'
%   • 'statTablehdrs'
%   • 'statTablerows': sr x 1 cell where each cell contains the label
%       associated with the data in the rows of statStruct
%   • 'statStruct': nel x 1 cell where each cell contains a 
%       sr x ntypeblks x (1+nregions) matrix where the rows correspond to
%       the values indicated in statTablerows, the columns correspond to
%       the statTablehdrs data and 3-rd dim corresonds to (1) all data, (2)
%       laser off @ beginning, (3) laser off @ end, (4) background or area
%       map, (5) tissue signal
%   • 'muperType':    ntypeblks x 4 x nel matrix with the mean value of
%       each element where the rows are ordered according to the appearance
%       of the type in the data structure d. The rows correspond to the
%       order of the labels in labelInfo.typerowInfo(:,4) (with names the
%       same as those in line_types(label)). Columns are organized as 
%   • 'stdperType'
%   • 'muperLine'
%   • 'medperLine'
%   • 'stdperLine'
%   • 'rmoutliers'
%   • 'lowerCO'
%   • 'upperCO'
%
%  Date           Author              E-mail                      Version
%   9 Oct  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1


if ~exist('hdrtxt')
  hdrtxt = [];
end

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

if strmatch('savestats',PropertyNames)
  savestats = PropertyVal{strmatch('savestats',PropertyNames)};
else
  savestats = 0;
end

if strmatch('rmout',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  rmoutliers = PropertyVal{strmatch('rmout',PropertyNames)};
else
  rmoutliers = 1;
end

if strmatch('normbyel',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  normbyel = PropertyVal{strmatch('normbyel',PropertyNames)};
else
  normbyel = [];
end

%---> Initialize workspace
delcell = struct2cell(rmfield(d,skip_fields));
colstart = 1; colend = 2; colLen = 3; colLabel = 4;
normbyeli = interp_data_elemrng(d,normbyel);
nel = length(labelInfo.elnames);
M = size(labelInfo.dchMask,1);
N = size(labelInfo.dchMask,2);
ncol = size(labelInfo.dchMask,3);
ntypeblks = size(labelInfo.typerowInfo,1);
ntypes = length(labelInfo.line_types);
stathdri = 0;
ztag = {'Whole Area';'Laser Off Beginning';'Laser Off End';'Background';'Tissue'};
statTablerows = {'Time';'# points';'min';'max';'mean';'median';'stnd. dev'};
statTablehdrs = cell(1,ntypeblks);
statStruct = cell(nel,1);
statStruct(:) = {zeros(10,ntypeblks,5)};
startdnum = datenum(d.line_datevec(1,:));
medperType = nan(ntypeblks,ncol,nel);
muperType = nan(ntypeblks,ncol,nel); % rows are types, columns are 1: LAbeginning, 2: LA end, 3: bkgrnd, 4: Tissue
stdperType = nan(ntypeblks,ncol,nel);
lowerCO =  nan(ntypes,ncol+1,nel);
upperCO =  nan(ntypes,ncol+1,nel);
muperLine = nan(M,ncol,nel);
medperLine = nan(M,ncol,nel);
stdperLine = nan(M,ncol,nel);

for t = 1:ntypes
  %---> Identify relevant parameters in datasetgiven type label
  blkrows = find(labelInfo.typerowInfo(:,colLabel) == t);
  nblks = length(blkrows);
  blkhr = (datenum(d.line_datevec(labelInfo.typerowInfo(blkrows,colend),:)) - startdnum).*24;
  typelabel = labelInfo.line_types{t};

  %---> Determine upper and lower cutoff values for outlier removal
  eltypeallcell = cellfun(@(x) x(labelInfo.type_rows{t},:),delcell,'uniformoutput',0);
  if ~isempty(normbyel)
    eltypeallcell = cellfun(@(x) x./delcell{normbyeli}(labelInfo.type_rows{t},:),eltypeallcell,'uniformoutput',0);
  end
  
  if rmoutliers
    for labelval = 0:4 
      if labelval == 0
        eltypeallcat = cellfun(@(x) x(labelInfo.dLabelMask(labelInfo.type_rows{t},:) > labelval),eltypeallcell,'uniformoutput',0);
      else
        eltypeallcat = cellfun(@(x) x(labelInfo.dLabelMask(labelInfo.type_rows{t},:) == labelval),eltypeallcell,'uniformoutput',0);
      end
      rawmean = cellfun(@(x) mean(x),eltypeallcat);
      allcatq1 = cellfun(@(x) quantile(x,0.25),eltypeallcat);
      allcatq3 = cellfun(@(x) quantile(x,0.75),eltypeallcat);
      allcatiqr = allcatq3 - allcatq1;
      allcatiqr(allcatiqr == 0) = rawmean(allcatiqr == 0)./3;
      lowerCO(t,labelval+1,:) =  allcatq1 - 3.*allcatiqr;
      upperCO(t,labelval+1,:) =  allcatq3 + 3.*allcatiqr;
    end
  else
    lowerCO(t,:,:) = repmat(cellfun(@(x) min(x(:)),eltypeallcell),1,5);
    upperCO(t,:,:) = repmat(cellfun(@(x) max(x(:)),eltypeallcell),1,5);
  end
  
  %---> Compute statistics and relevant features in each block for each
  %type 
  for b = 1:nblks
    stathdri = stathdri + 1;
    if ~isempty(strfind(lower(typelabel),'calib')) || ~isempty(strfind(lower(typelabel),'nist'))
      blkr = labelInfo.typerowInfo(blkrows(b),1):labelInfo.typerowInfo(blkrows(b),2);
      blkmask = labelInfo.dLabelMask(blkr,:);%tmask((sum(labelInfo.typerowInfo(blkrows(1:b-1),3))+1):sum(labelInfo.typerowInfo(blkrows(1:b),3)),:);
      statTablehdrs{stathdri} = strcat(typelabel,'-',num2str(b));
    else
      blkr = labelInfo.type_rows{t};
      blkmask = labelInfo.dLabelMask(blkr,:);%tmask((sum(labelInfo.typerowInfo(blkrows(1:b-1),3))+1):sum(labelInfo.typerowInfo(blkrows(1:b),3)),:);
      statTablehdrs{stathdri} = strcat(typelabel,'-',num2str(length(blkrows)));
    end
    
    for el = 1:nel
      for cdat = 1:5
        labelval = cdat - 1;
        eldat = delcell{el}(blkr,:);
        if ~isempty(normbyel)
          eldat = eldat./delcell{normbyeli}(blkr,:);
        end
      
        if cdat > 1
          nanmat = nan(size(blkmask));
          nanmat(blkmask == cdat - 1) = 1;
          nblkpnts = sum(sum(blkmask == cdat - 1));
          nanmat(eldat < lowerCO(t,cdat,el)) = nan;
          if sum(eldat(blkmask == cdat - 1) > upperCO(t,cdat,el)) < nblkpnts*.5 % most points are less than the upper cutoff
            nanmat(eldat > upperCO(t,cdat,el)) = nan;
          end
          celdat = eldat.*nanmat;
          muperLine(blkr,labelval,el) = nanmean(celdat,2);
          medperLine(blkr,labelval,el) = nanmedian(celdat,2);
          stdperLine(blkr,labelval,el) = nanstd(celdat,[],2);
          
          if ~isempty(strfind(lower(typelabel),'calib')) || ~isempty(strfind(lower(typelabel),'nist'))
            medperType(blkrows(b),labelval,el) = nanmedian(celdat(:));
            muperType(blkrows(b),labelval,el) = nanmean(celdat(:));
            stdperType(blkrows(b),labelval,el) = nanstd(celdat(:));
          else
            for sb = 1:length(blkrows)
              sbblkr = labelInfo.typerowInfo(blkrows(sb),1):labelInfo.typerowInfo(blkrows(sb),2);
              sbblkmask = labelInfo.dLabelMask(sbblkr,:);
              elsb = delcell{el}(sbblkr,:);
              if ~isempty(normbyel)
                elsb = elsb./delcell{normbyeli}(sbblkr,:);
              end
              nanmat = nan(size(sbblkmask));
              nanmat(sbblkmask == cdat - 1) = 1;
              nanmat(elsb < lowerCO(t,cdat,el)) = nan;
              nanmat(elsb > upperCO(t,cdat,el)) = nan;
              celdat = elsb.*nanmat;
              muperType(blkrows(sb),labelval,el) = nanmean(celdat(:));
              stdperType(blkrows(sb),labelval,el) = nanstd(celdat(:));
            end
          end
        end
        
        switch labelval
          case 0
            eldat = eldat(blkmask > labelval);
          otherwise
            eldat = eldat(blkmask == labelval);
        end
         eldat = eldat(eldat >= lowerCO(t,cdat,el) & eldat <= upperCO(t,cdat,el));
        
        if ~isempty(eldat)
          coli = 1;
          statStruct{el}(coli,stathdri,cdat) = blkhr(b); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = length(eldat(:)); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = min(eldat); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = max(eldat); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = mean(eldat); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = median(eldat); coli = coli + 1;
          statStruct{el}(coli,stathdri,cdat) = std(eldat); coli = coli + 1;
        end %---> end stat loop
      end %---> end of class type loop (ALL, LA off before, LA off end, background, tissue)
    end %---> end element loop
  end %---> end block loop
end %---> end type loop

dstats.ntypeblks = ntypeblks;
dstats.statTablehdrs = statTablehdrs;
dstats.statTablerows = statTablerows;
dstats.statStruct = statStruct;
dstats.medperType = medperType;
dstats.muperType = muperType;
dstats.stdperType = stdperType;
dstats.medperLine = medperLine;
dstats.muperLine = muperLine;
dstats.stdperLine = stdperLine;
dstats.rmoutliers = rmoutliers;
dstats.lowerCO = lowerCO;
dstats.upperCO = upperCO;

if savestats
  cdir = loaddirfun;
  xlpath = strcat(cdir.fpath,'MIMS Analysis\',d.dataset,' Region Stats - 31{ ratio, outlier adjusted.xlsx');
  ntabc = sum(cellfun(@(x) ~isempty(x),statTablehdrs)) + 1;
  ntabr = length(statTablerows)+1;
  for el = 1:nel
    elstr = elnames{el};
    elstat = statStruct{el};
    dcell = cell((ntabr+2) * size(elstat,3),ntabc);
    for c = 1:size(elstat,3)
      fillrow = (c-1)*(ntabr+2) + 1;
      dcell{fillrow,1} = ztag{c};
      dcell((fillrow+1):(fillrow + ntabr - 1),1) = statTablerows;
      dcell(fillrow,2:ntabc) = statTablehdrs(1:ntabc-1);
      dcell((fillrow+1):(fillrow + ntabr - 1),2:ntabc) =  num2cell(statStruct{el}(1:ntabr-1,1:ntabc-1,c));
    end
    xlswrite(xlpath,dcell,elstr);
  end
  sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
  Excel = actxserver('Excel.Application');
  Excel.Workbooks.Open(fullfile(xlpath)); % Full path is necessary!
  try
    Excel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
  catch
    ; % Do nothing.
  end
  % Save, close and clean up.
  Excel.ActiveWorkbook.Save;
  Excel.ActiveWorkbook.Close;
  Excel.Quit;
  Excel.delete;
end

end