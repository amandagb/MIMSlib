function data_struct = readIVISfiles(varargin);
%% Script:  data_struct = readIVISfiles(varargin);
% d = readIVISfiles('dat','20141107_nc4V_ivisbs');
% Description:  Reads files from the IVIS software. Parses descriptive text
% files and collects information about each of the images.
% Example:  d = readfiles;
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   • 'data_dir':  string indicating the data directory [DEFAULT = user selected]
%                  EXAMPLE: 'MAD20141107150558_SEQ'
%   • 'clicknum': vector indicating which click numbers to load [DEFAULT =
%   all]
%   • 'unmixed': logical indicating whether to load the spectral unmixed
%         data only (1) or not (0) [DEFAULT = 0];
%   • 'photo': logical indicating whether to load the spectral unmixed
%         data only (1) or not (0) [DEFAULT = 1];
%   • 'lum': logical indicating whether to load the spectral unmixed
%         data only (1) or not (0) [DEFAULT = 1];
%   • 'fluorref': logical indicating whether to load the spectral unmixed
%         data only (1) or not (0) [DEFAULT = 1];
%   • 'readbias': logical indicating whether to load the spectral unmixed
%         data only (1) or not (0) [DEFAULT = 1];
%   • 'allphotos': logical indicating whether to include all photos from
%         the individual folders [DEFAULT = 0];
% OUTPUTS ---------------------------------------------------------------
% data_struct:  a structure contining a field with the name of each
%   element for which data was collectedata_struct. Each field contains a <# lines> x
%   <# time samples> double matrix.
%  Date           Author              E-mail                      Version
%  18 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%  19 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1
%     createtable function lists fileds in alphabetical order
%  23 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1.1
%     changed createtable to list fields in order of occurence in the text
%     file
%  24 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1.2
%     Changing the way the images are contained
%  26 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     2
%     Reflects new naming scheme -- does not require the parent data
%     directory to influences the call to childen directories
%  23 Feb  2015   Amanda Balderrama   amanda.gaudreau@gmail.com     3
%     Updated to handle new naming scheme introduced by adding spectral
%     unmixed data to folders
%  27 Feb  2015   Amanda Balderrama   amanda.gaudreau@gmail.com     3.1
%  9  Mar  2015   Amanda Balderrama   amanda.gaudreau@gmail.com     3.2
%  12 Mar  2015   Amanda Balderrama   amanda.gaudreau@gmail.com     3.3
%     Reads in now "FloatCorrected" versions of tiffs if available

%% ---- Read data from text files, assumes comma delimiter ----
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

cpaths = loaddirfun;
current_dir = cpaths.ivisfldr; %MADLab Research\Data\

if strmatch('dat',PropertyNames)
  v = PropertyVal{strmatch('dat',PropertyNames)};
  dir_name = strcat(current_dir,v);
else
  dir_name = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end
while ~isdir(dir_name)
  dir_name = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 1;
end

if isempty(strmatch('save_data',PropertyNames)) & strmatch('savedata',PropertyNames)
  svdat = PropertyVal{strmatch('savedata',PropertyNames)};
end

if strmatch('unmixed',PropertyNames)
  umxonly = PropertyVal{strmatch('unmixed',PropertyNames)};
else
  umxonly = 0;
end

if strmatch('photo',PropertyNames)
  incphoto = PropertyVal{strmatch('photo',PropertyNames)};
else
  incphoto = 1;
end

if strmatch('allphoto',PropertyNames)
  allphoto = PropertyVal{strmatch('allphoto',PropertyNames)};
else
  allphoto = 0;
end

if strmatch('lum',PropertyNames)
  inclum = PropertyVal{strmatch('lum',PropertyNames)};
else
  inclum = 1;
end

if strmatch('readbias',PropertyNames)
  increadbias = PropertyVal{strmatch('readbias',PropertyNames)};
else
  increadbias = 1;
end

if strmatch('fluorref',PropertyNames)
  incfluorref = PropertyVal{strmatch('fluorref',PropertyNames)};
else
  incfluorref = 1;
end

if strmatch('clicknum',PropertyNames)
  clkvec = PropertyVal{strmatch('clicknum',PropertyNames)};
else
  clkvec = [];
end

cd (dir_name); % got to that directory
slashi = strfind(dir_name,'\');
data_fldr = dir_name(slashi(end)+1:end);
usi = strfind(data_fldr,'_'); % underscore indices
data_header = data_fldr(1:usi-1);
allfiles = ls;
dirlistcell = cellstr(allfiles(3:end,:));
[hdrs,destxt] = strtok(dirlistcell,'_');
[fileorfldr,ftype] = strtok(dirlistcell,'.');
folderbin = cellfun(@(x) isempty(x),ftype);
fldrlist = dirlistcell(folderbin);
[~,fldrending] = strtok(fldrlist,'_');
fldrending = strrep(fldrending,'_','');
seqfile = strfind(dirlistcell,strcat('SequenceInfo.txt'));
seqf_indx = find(cellfun(@(x) ~isempty(x) && x == 1,seqfile));
sequmxf_indx = find(cellfun(@(x) ~isempty(x) && x > 1,seqfile));

if ~isempty(seqf_indx)
  seqfid = fopen(dirlistcell{seqf_indx});
  [seqhtxt,seqftxt] = readtxtlines(seqfid);
  fclose(seqfid);
  data_struct.SequenceInfo.headers = cutcell(seqhtxt);
  data_struct.SequenceInfo.info = cutcell(seqftxt);
  clkind = nonemptystrind(data_struct.SequenceInfo.info(:,2),'ClickNumber ');
  nclks = length(clkind);
else
  disp('No SequenceInfo.txt file available for this sequence');
  nclks = length(fldrlist);
end

if ~isempty(sequmxf_indx) || isempty(seqf_indx)
  if length(sequmxf_indx) > 1
    for seqi = 1:length(sequmxf_indx)
      sequmxfid = fopen(dirlistcell{sequmxf_indx(seqi)});
      [sequmxhtxt,sequmxftxt] = readtxtlines(sequmxfid);
      fclose(sequmxfid);
      data_struct.SequenceInfoUMX.headers{seqi} = cutcell(sequmxhtxt);
      data_struct.SequenceInfoUMX.info{seqi} = cutcell(sequmxftxt);
    end
  else
    sequmxfid = fopen(dirlistcell{sequmxf_indx});
    [sequmxhtxt,sequmxftxt] = readtxtlines(sequmxfid);
    fclose(sequmxfid);
    data_struct.SequenceInfoUMX.headers = cutcell(sequmxhtxt);
    data_struct.SequenceInfoUMX.info = cutcell(sequmxftxt);
    if isempty(seqf_indx)
      data_struct.SequenceInfo.headers = cutcell(sequmxhtxt);
      data_struct.SequenceInfo.info = cutcell(sequmxftxt);
    end
  end
end

umxfldrind = find(cellfun(@(x) ~isempty(x),strfind(fldrlist,'UMX')));
fldrnums = cellfun(@(x) str2num(x),fldrending,'uniformoutput',0);
numfldrind = find(cellfun(@(x) ~isempty(x),fldrnums));
fldrnums = cell2mat(fldrnums(numfldrind));
foffset = max(numfldrind) - length(fldrnums);

if umxonly
  usefldrbin = umxfldrind;
elseif ~isempty(clkvec)
  usefldrbin = numfldrind(ismember(fldrnums,clkvec)); 
  if any(clkvec > length(fldrlist))
    usefldrbin = [];
  elseif any(clkvec > max(fldrnums))
    usefldrbin = [usefldrbin(:)',clkvec(clkvec > max(fldrnums)) + foffset];
  end
else
  [Ustr,ic,ia] = unique(strtok(fldrlist,'_'));
  lstru = zeros(1,length(Ustr));
  for u = 1:length(Ustr)
    lstru(u) = sum(u == ia);
  end
  [~,usestri] = max(lstru);
  usefldrbin = find(ia == usestri);
end

% Open subfolders within data_header
imfields = {'fluorescent reference','luminescent','luminescentFloatCorrected',...
  'photograph','photographFloatCorrected','readbiasonly'};
imfldsnospace = strrep(imfields,' ','');
fstr = 'fluorescentreference'; 
frfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);
fstr = 'luminescent'; 
lumfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);
fstr = 'luminescentFloatCorrected'; 
lumcorrfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);
fstr = 'photograph'; 
pfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);
fstr = 'photographFloatCorrected'; 
pcorrfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);
fstr = 'readbiasonly'; 
rbfillcol = cellfun(@(x,y) ~isempty(x) && length(y) == length(fstr),strfind(imfldsnospace,fstr),imfldsnospace);


data_struct.imageNames = cell(length(usefldrbin)+1,length(imfields));
data_struct.imageNames(1,:) = imfields;
if isempty(usefldrbin)
  disp(sprintf('No folders match user qeury for data folder input %s',data_fldr));
else
  for fs = 1:length(usefldrbin)
    cd(fldrlist{usefldrbin(fs)}) %cd(data_struct.SequenceInfo.info{clkind(fs),3})
    disp(sprintf('Processing file %s',fldrlist{usefldrbin(fs)}))
    
    % Analyzed Click Info Text File
    anclkfid = fopen('AnalyzedClickInfo.txt');
    [anclkihtxt,anclkiftxt] = readtxtlines(anclkfid);
    fclose(anclkfid);
    data_struct.AnalyzedClickInfo.headers{fs} = cutcell(anclkihtxt);
    data_struct.AnalyzedClickInfo.info{fs} = cutcell(anclkiftxt);
    if fs == 1
      data_struct.AnalyzedTable = data_struct.AnalyzedClickInfo.info{fs};
    else
      data_struct.AnalyzedTable = createtable(data_struct.AnalyzedTable,...
        data_struct.AnalyzedClickInfo.info{fs});
    end
    
    % Click Info Text File
    clkfid = fopen('ClickInfo.txt');
    [clkihtxt,clkiftxt] = readtxtlines(clkfid);
    fclose(clkfid);
    data_struct.ClickInfo.headers{fs} = cutcell(clkihtxt);
    data_struct.ClickInfo.info{fs} = cutcell(clkiftxt);
    if fs == 1
      data_struct.ClickTable = data_struct.ClickInfo.info{fs};
    else
      data_struct.ClickTable = createtable(data_struct.ClickTable,...
        data_struct.ClickInfo.info{fs});
    end
    
    if incfluorref
      fstr = 'fluorescentreference';
      if exist(strcat(fstr,'.tif')) == 2
        fluorescentreference.image{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,frfillcol} = fldrlist{usefldrbin(fs)};
      else
        fluorescentreference.image{fs} = [];
      end
    else
      fluorescentreference.image{fs} = [];
    end
    
    if inclum
      fstr = 'luminescent';
      if exist(strcat(fstr,'.tif')) == 2
        luminescent.image{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,lumfillcol} = fldrlist{usefldrbin(fs)};
      else
        luminescent.image{fs} = [];
      end
    else
      luminescent.image{fs} = [];
    end
    
    if inclum
      fstr = 'luminescentFloatCorrected';
      if exist(strcat(fstr,'.tif')) == 2
        luminescent.imageCorrected{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,lumcorrfillcol} = fldrlist{usefldrbin(fs)};
      else
        luminescent.imageCorrected{fs} = [];
      end
    else
      luminescent.imageCorrected{fs} = [];
    end
    
    if incphoto
      fstr = 'photograph';
      if exist(strcat(fstr,'.tif')) == 2 && allphoto == 1
        photograph.image{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,pfillcol} = fldrlist{usefldrbin(fs)};
      elseif exist(strcat(fstr,'.tif')) == 2 && fs == 1
        photograph.image{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,pfillcol} = fldrlist{usefldrbin(fs)};
      else
        photograph.image{fs} = [];
      end
    else
      photograph.image{fs} = [];
    end
    
    if incphoto
      fstr = 'photographFloatCorrected';
      if exist(strcat(fstr,'.tif')) == 2 && allphoto == 1
        photograph.imageCorrected{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,pcorrfillcol} = fldrlist{usefldrbin(fs)};
      elseif exist(strcat(fstr,'.tif')) == 2 && fs == 1
        photograph.imageCorrected{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,pcorrfillcol} = fldrlist{usefldrbin(fs)};
      else
        photograph.imageCorrected{fs} = [];
      end
    else
      photograph.imageCorrected{fs} = [];
    end
    
    if increadbias
      fstr = 'readbiasonly';
      if exist(strcat(fstr,'.tif')) == 2
        readbiasonly.imageCorrected{fs} = imread(strcat(fstr,'.tif'));
        data_struct.imageNames{fs + 1,rbfillcol} = fldrlist{usefldrbin(fs)};
      else
        readbiasonly.image{fs} = [];
      end
    else
      readbiasonly.image{fs} = [];
    end
    cd(dir_name)
  end
  
  data_struct.ClickTable = cutcell(data_struct.ClickTable);
  data_struct.AnalyzedTable = cutcell(data_struct.AnalyzedTable);
  
  
  if any(cellfun(@(x) ~isempty(x),fluorescentreference.image))
    fluorescentreference = createIstruct(fluorescentreference,imfields{frfillcol},...
      data_struct.ClickInfo.headers,data_struct.ClickTable,...
      data_struct.AnalyzedClickInfo.headers,data_struct.AnalyzedTable);
    data_struct.fluorescentreference = fluorescentreference;
  end
  
  if any(cellfun(@(x) ~isempty(x),luminescent.image))
    luminescent = createIstruct(luminescent,imfields{lumfillcol},...
      data_struct.ClickInfo.headers,data_struct.ClickTable,...
      data_struct.AnalyzedClickInfo.headers,data_struct.AnalyzedTable);
    data_struct.luminescent = luminescent;
  end
  
  if any(cellfun(@(x) ~isempty(x),photograph.image))
    photograph = createIstruct(photograph,imfields{pfillcol},...
      data_struct.ClickInfo.headers,data_struct.ClickTable,...
      data_struct.AnalyzedClickInfo.headers,data_struct.AnalyzedTable);
    data_struct.photograph = photograph;
  end
  
  if any(cellfun(@(x) ~isempty(x),readbiasonly.image))
    readbiasonly = createIstruct(readbiasonly,imfields{rbfillcol},...
      data_struct.ClickInfo.headers,data_struct.ClickTable,...
      data_struct.AnalyzedClickInfo.headers,data_struct.AnalyzedTable);
    data_struct.readbiasonly = readbiasonly;
  end
  
  seqind = nonemptystrind(data_struct.SequenceInfo.info(:,2),'Sequence');
  data_struct.origseq = data_struct.SequenceInfo.info{seqind,3};
  data_struct.dataset = data_fldr;
end



%=====================================================
% SUBFUNCTION
%=====================================================
  function [hdrtxt,fldtxt] = readtxtlines(fid)
    linecnt = 0;
    header_lines = [];
    hdrtxt = cell(60,3); ht = 1;
    fldtxt = cell(1200,4); ft = 1;startch = 0;
    tline = fgetl(fid);
    while ischar(tline)
      linecnt = linecnt + 1;
      colons = strfind(tline,':');
      hstars = strfind(tline, '***');
      cmnt = strfind(tline,'#');
      hcmnt = strfind(tline,'#################################################################');
      if hstars
        hdrtxt{ht,1} = linecnt;
        hdrtxt{ht,2} = tline(5:colons-1);
        hdrtxt{ht,3} = 's';
        fldtxt{ft,1} = linecnt;
        fldtxt{ft,2} = strcat('*>',hdrtxt{ht,2});
        fldtxt{ft,3} = tline(colons+2:end);
        ht = ht + 1;
        ft = ft + 1;
      elseif hcmnt
        if startch
          startch = 0;
        else
          startch = 1;
        end
      elseif cmnt
        if startch
          hdrtxt{ht,1} = linecnt;
          hdrtxt{ht,2} = tline(cmnt+1:end);
          hdrtxt{ht,3} = 'p';
          fldtxt{ft,1} = linecnt;
          fldtxt{ft,2} = strcat('#>',hdrtxt{ht,2});
          ft = ft + 1;
          ht = ht + 1;
        elseif any(cmnt == 1)
          fldtxt{ft,1} = linecnt;
          fldtxt{ft,2} = tline;
          ft = ft + 1;
        elseif isempty(colons)
          fldtxt{ft,1} = linecnt;
          fldtxt{ft,2} = tline;
          ft = ft + 1;
        elseif length(colons) > 0
          fldtxt{ft,1} = linecnt;
          fldtxt{ft,2} = sprintf('%s %s',tline(1:colons(1)-1),tline(cmnt(1):end));
          fldtxt{ft,3} = tline(colons(1)+1:cmnt-1);
          ft = ft + 1;
        end
      elseif colons
        fldtxt{ft,1} = linecnt;
        fldtxt{ft,2} = tline(1:colons-1);
        fldtxt{ft,3} = tline(colons+2:end);
        ft = ft + 1;
      end
      tline = fgetl(fid);
    end
  end

%=====================================================
% SUBFUNCTION
%=====================================================
  function cellnoemptyrows = cutcell(cellwempty);
    nonlines = cellfun(@(x) ~isempty(x),cellwempty);
    keep_rows = find(sum(nonlines(:,2:end),2) ~= 0);
    cellnoemptyrows = cellwempty(keep_rows,:);
    keep_cols = find(sum(nonlines,1) ~= 0);
    cellnoemptyrows = cellnoemptyrows(:,keep_cols);
  end

%=====================================================
% SUBFUNCTION
%=====================================================
  function outcell = createtable(currentcell,newcell)
    currentcell = cutcell(currentcell);
    newcell = cutcell(newcell);
    [cm,cn] = size(currentcell);
    currentstr = cellstr(currentcell(:,2));
    newstr = cellstr(newcell(:,2));
    newhdrind = sort([find(cellfun(@(x) ~isempty(x),strfind(newstr,'*>')));...
      find(cellfun(@(x) ~isempty(x),strfind(newstr,'#>')))]);
    newhdr = newstr(newhdrind);
    
    currhdrind = sort([find(cellfun(@(x) ~isempty(x),strfind(currentstr,'*>')));...
      find(cellfun(@(x) ~isempty(x),strfind(currentstr,'#>')))]);
    currhdr = currentstr(currhdrind);
    
    [~,indcs,indns] = intersect(currhdr,newhdr,'stable');
    [~,ixorcs,ixorns] = setxor(currhdr,newhdr,'stable');
    
    lenc = length(currentstr);
    lenn = length(newstr);
    totalrows = min(lenn,lenc);
    outcell = cell(totalrows,size(currentcell,2)+1);
    
    lasti = 1;
    for c = 1:length(indcs)
      newind = indns(c);
      newstarti = newhdrind(newind);
      if newstarti == lenn || newind == length(newhdrind)
        newendi = lenn;
      else
        newendi = newhdrind(newind+1)-1;
      end
      newirng = newstarti:newendi;
      
      currind = indcs(c);
      currstarti = currhdrind(currind);
      if currstarti == lenc || currind == length(currhdrind)
        currendi = lenc;
      else
        currendi = currhdrind(currind+1)-1;
      end
      currirng = currstarti:currendi;
      
      newhdrfields = newcell(newirng,2);
      currhdrfields = currentcell(currirng,2);
      
      [commonflds,indfcs,indfns] = intersect(currhdrfields,newhdrfields,'stable');
      [diffflds,dics,dins] = setxor(currhdrfields,newhdrfields,'stable');
      
      % --- Fill existing fields
      outrows = lasti:(length(indfcs)+lasti-1);
      outcell(outrows,1) = ...
        mat2cell([[cell2mat(currentcell(currirng(indfcs),1))],...
        [cell2mat(newcell(newirng(indfns),1))]],...
        [ones(1,length(indfcs))]);
      outcell(outrows,2) = commonflds;
      outcell(outrows,3:cn) = currentcell(currirng(indfcs),3:cn);
      outcell(outrows,end) = newcell(newirng(indfns),end);
      lasti = outrows(end) + 1;
      
      if ~isempty(dics)
        outrows = lasti:(length(dics)+lasti-1);
        outcell(outrows,1) = ...
          mat2cell([[cell2mat(currentcell(currirng(dics),1))],...
          [zeros(length(outrows),1)]],...
          [ones(1,length(dics))]);
        outcell(outrows,2:cn) = currentcell(currirng(dics),2:cn);
        lasti = outrows(end) + 1;
      end
      
      if ~isempty(dins)
        outrows = lasti:(length(dins)+lasti-1);
        outcell(outrows,1) = ...
          mat2cell([[zeros(length(dins),cn-2)],...
          [cell2mat(newcell(newirng(dins),1))]], [ones(1,length(dins))]);
        outcell(outrows,2) = newhdrfields(dins);
        outcell(outrows,end) = newcell(newirng(dins),end);
        lasti = outrows(end) + 1;
      end
    end
    
    if ~isempty(ixorcs)
      for c = 1:length(ixorcs)
        currind = ixorcs(c);
        currstarti = currhdrind(currind);
        if currstarti == lenc || currind == length(currhdrind)
          currendi = lenc;
        else
          currendi = currhdrind(currind+1)-1;
        end
        currirng = currstarti:currendi;
        
        outrows = lasti:(length(currirng)+lasti-1);
        outcell(outrows,1) = ...
          mat2cell([[cell2mat(currentcell(currirng,1))],...
          [zeros(length(outrows),1)]],...
          [ones(1,length(outrows))]);
        outcell(outrows,2:cn) = currentcell(currirng,2:cn);
        lasti = outrows(end) + 1;
      end
    end
    
    if ~isempty(ixorns)
      for c = 1:length(ixorns)
        newind = ixorns(c);
        newstarti = newhdrind(newind);
        if newstarti == lenn || newind == length(newhdrind)
          newendi = lenn;
        else
          newendi = newhdrind(newind+1)-1;
        end
        newirng = newstarti:newendi;
        
        outrows = lasti:(length(newirng)+lasti-1);
        outcell(outrows,1) = ...
          mat2cell([[zeros(length(outrows),cn-2)],...
          [cell2mat(newcell(newirng,1))]], [ones(1,length(outrows))]);
        outcell(outrows,[2,end]) = newcell(newirng,2:end);
        lasti = outrows(end) + 1;
      end
    end
    
    outcell = cutcell(outcell);
  end

%=====================================================
% SUBFUNCTION
%=====================================================
  function Istruct = createIstruct(Istruct,imgstr,headersct,tablect,headersan,tablean)
    [~,hmaxelct] = max(cellfun(@(x) size(x,1),headersct));
    
    hstrct = strfind(headersct{hmaxelct}(:,2),imgstr);
    hindct = find(cellfun(@(x) ~isempty(x),hstrct));
    fullhstrct = headersct{hmaxelct}{hindct,2};
    
    hstr_tabct = strfind(tablect(:,2),fullhstrct);
    hind_tabct = find(cellfun(@(x) ~isempty(x),hstr_tabct))+1;
    
    nexthstr_tabct = strfind(tablect(:,2),'>');
    nexthv_tabct = find(cellfun(@(x) ~isempty(x),nexthstr_tabct));
    nexthvposct = find(nexthv_tabct > hind_tabct);
    nexthind_tabct = nexthv_tabct(nexthvposct(1))-1;
    
    [~,hmaxelan] = max(cellfun(@(x) size(x,1),headersan));
    
    hstran = strfind(headersan{hmaxelan}(:,2),imgstr);
    hindan = find(cellfun(@(x) ~isempty(x),hstran));
    fullhstran = headersan{hmaxelan}{hindan,2};
    
    hstr_taban = strfind(tablean(:,2),fullhstran);
    hind_taban = find(cellfun(@(x) ~isempty(x),hstr_taban))+1;
    
    nexthstr_taban = strfind(tablean(:,2),'>');
    nexthv_taban = find(cellfun(@(x) ~isempty(x),nexthstr_taban));
    nexthvposan = find(nexthv_taban > hind_taban);
    nexthind_taban = nexthv_taban(nexthvposan(1))-1;
    
    imgheadersct = tablect(hind_tabct:nexthind_tabct,2);
    imgheadersan = tablean(hind_taban:nexthind_taban,2);
    if length(imgheadersct) ~= length(imgheadersan)
      [minL,CTorAN] = min([length(imgheadersct),length(imgheadersan)]);
      if any(strcmp(imgheadersct(1:minL),imgheadersan(1:minL)) == 0)
        disp('Difference between Click and AnalyzedClick information files');
      elseif CTorAN == 1
        imgheaders = imgheadersan;
        hind_tab = hind_taban;
        nexthind_tab = nexthind_taban;
        table = tablean;
      else
        imgheaders = imgheadersct;
        hind_tab = hind_tabct;
        nexthind_tab = nexthind_tabct;
        table = tablect;
      end
    elseif any(strcmp(imgheadersct,imgheadersan) == 0)
      disp('Difference between Click and AnalyzedClick information files');
    else
      imgheaders = imgheadersct;
      hind_tab = hind_tabct;
      nexthind_tab = nexthind_tabct;
      table = tablect;
    end
    
    for ih = hind_tab:nexthind_tab
      dline = table(ih,:);
      nel = length(dline)-2;
      nonemptyd = find(cellfun(@(x) ~isempty(x),dline(3:end)))+2;
      if ~isempty(nonemptyd)
        dcontent = str2num(dline{nonemptyd(end)});
        nnonempt = length(nonemptyd);
        if all(cellfun(@(x) isempty(x), dline(3:end)))
          Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),[]);
        elseif isempty(dcontent)
          fillcell = cell(1,nel);
          fillcell(nonemptyd-2) = dline(nonemptyd);
          Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),fillcell);
        elseif ~isempty(strmatch(dline{nonemptyd(end)},num2str(dcontent)))
          fillcell = nan(1,nel);
          fillcell(nonemptyd-2) = str2num(char(dline(nonemptyd)))';
          Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),fillcell);
        elseif isempty(dline{3})
          fillcell = cell(1,nel);
          fillcell(nonemptyd-2) = dline(nonemptyd);
          Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),fillcell);
        else
          fillcell = cell(1,nel);
          fillcell(nonemptyd-2) = dline(nonemptyd);
          Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),fillcell);
        end
      else
        Istruct = setfield(Istruct,regexprep(dline{2},'\W',''),cell(1,nel));
      end
    end
    
    nimgs = length(Istruct.image);
    cnt2phC = zeros(1,nimgs);
    for i = 1:length(imgstr)
      Ifields = fieldnames(Istruct);
      FOVid = nonemptystrind(lower(Ifields),'fieldofview');
      if ~isempty(FOVid)
        FOVnum = getfield(Istruct,Ifields{FOVid});
        fnumid = nonemptystrind(lower(Ifields),'fnumber');
        fnumnum = getfield(Istruct,Ifields{fnumid});
        
        for n = 1:nimgs
          if rem(FOVnum(n),1) ~= 0
            cnt2phCind = nonemptystrind(tablect(:,2),sprintf('Coef C-ccd at FOV %1.1f, f%d',FOVnum(n),fnumnum(n)));
          else
            cnt2phCind = nonemptystrind(tablect(:,2),sprintf('Coef C-ccd at FOV %d, f%d',FOVnum(n),fnumnum(n)));
          end
          cnt2phC(n) = str2num(tablect{cnt2phCind,3});
        end
      end
      Istruct = setfield(Istruct,'cnts2photonC',cnt2phC);
    end
  end
end
