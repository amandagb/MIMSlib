function [varargout] = parseCBMstruct(d,varargin);
%% Script:  [varargout] = parseCBMstruct(d,varargin);
% Description:
% Example:
% Required Function: loaddirfun, skip_fields, auto_thresh, save32bitTIF
% INPUTS ----------------------------------------------------------------
% data_struct:  a structure contining a field with the name of each
%   element for which data was collected. Each field contains a <# lines> x
%   <# time samples> double matrix.
% varargin - 'PropertyName','PropertyValue'
%   • 'nest':  logical indicating whether to nest the data within the
%   parent directory data folder (1) or not (0). If 0, then a new folder
%   within the laser ablation folder will be created with the standard file
%   structure [DEFAULT = 0; 1 for NC]
% OUTPUTS ---------------------------------------------------------------
% NONE (as of v0)
%  Date           Author            E-mail                      Version
%   4 Mar  2015   Amanda Balderrama amanda.gaudreau@gmail.com     0
%   6 Mar  2015   Amanda Balderrama amanda.gaudreau@gmail.com     0.1
%  29 Mar  2015   Amanda Balderrama amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

compdir = loaddirfun;
if strmatch('comp',PropertyNames)
  v = PropertyVal{strmatch('comp',PropertyNames)};
  if strmatch(v,'pc')
    current_dir = compdir.lafldr; % 'C:\Users\ADGB\Documents\MADLab Data\Laser Ablation\';
  elseif strmatch(v,'pho')
    current_dir = 'C:\Users\MADLAB\Desktop\Laser Ablation\';
  elseif ~isempty(v)
    current_dir = v;
  else
    current_dir = cd;
  end
else
  current_dir = compdir.lafldr; % 'C:\Users\ADGB\Documents\MADLab Data\Laser Ablation\';
end

if ~isdir(current_dir)
  current_dir = cd;
end

if isempty(strfind(current_dir,'ADGB')) || isempty(strmatch('run4noel',PropertyNames))
  noelinit = 1;
else
  noelinit = 0;
end

if strmatch('nest',PropertyNames)
  nest = PropertyVal{strmatch('nest',PropertyNames)};
elseif noelinit
  nest = 0;
else
  nest = 0;
end

if strmatch('runread',PropertyNames)
  runread = PropertyVal{strmatch('runread',PropertyNames)};
elseif noelinit
  runread = 1;
else
  runread = 0;
end

if isstruct(d)
  origpath = strcat(current_dir,d.dataset);
  datasetfield = d.dataset;
elseif isstr(d)
  origpath =  strcat(current_dir,d);
  datasetfield = d;
else
  disp('Invalid input')
  return
end

cd(origpath);
fcontent = cellstr(ls);
fcontent = fcontent(3:end);
fldrind = cellfun(@(x) isempty(strfind(x,'.')),fcontent);
fileind = abs(fldrind - 1);

cd('RawData')
rawdatacontent = cellstr(ls);
rawdatacontent = rawdatacontent(3:end);
allhave = 0;
[stend,stri] = min(cellfun(@(x) length(x),rawdatacontent));
stend = stend - 4;
while ~allhave
  startstr = strfind(rawdatacontent,rawdatacontent{stri}(1:stend));
  strlog = cellfun(@(x) ~isempty(x),startstr);
  if sum(strlog) == length(rawdatacontent)
    allhave = 1;
  else
    stend = stend - 1;
  end
end
rawhdrstr = rawdatacontent{stri}(1:stend);
rawdatafmt = rawdatacontent{stri}(end-3:end);
linenum = cellfun(@(x) str2num(strrep(x(stend+1:end),rawdatafmt,'')),rawdatacontent);
[~,csvorder] = sort(linenum);

hdrfile = ismember(fcontent,'headerinfo.xlsx');
if any(hdrfile)
  [A,B,hdrtxt] = xlsread(strcat(origpath,'\',fcontent{hdrfile}));
else
  [~,hdrtxt] = readfiles('dat',datasetfield,'readcsv',0);
end

rownum = sum(cellfun(@(x) isnumeric(x),hdrtxt),2) >= 83;
allnanrow = find(sum(cellfun(@(x) isnan(x),hdrtxt(rownum,1:83)),2) >= 83);
rownumid = find(rownum);
rmrownum = rownumid(allnanrow);
keephdrrows = true(1,size(hdrtxt,1));
keephdrrows(rmrownum) = 0;
hdrtxt = hdrtxt(keephdrrows,:);

lineCO = strfind(lower(hdrtxt(2:end,1)),'_line');
relhtxt = cellfun(@(x,y) x(1:(y-1)),hdrtxt(2:end,1),lineCO,'uniformoutput',0);
reptxtind = find(cellfun(@(x) isempty(x),relhtxt));
relhtxt(reptxtind) = hdrtxt(reptxtind+1,1);
[line_types,tst,linesind] = unique(relhtxt,'stable');
[~,orderlines] = sort(tst);

if length(orderlines) > 1
  caliblines = find(cellfun(@isempty,strfind(lower(hdrtxt(2:end,1)),'calib')) == 0);
  sollines =  find(cellfun(@isempty,strfind(lower(hdrtxt(2:end,1)),'solution')) == 0);
  if ~isempty(caliblines)
    nieghborlabels = unique(relhtxt([max(caliblines(1)-1,1),min(caliblines(end)+1,length(relhtxt))],1));
    inccaliblines = find(ismember(lower(relhtxt),lower(nieghborlabels)));
    inccaliblines = inccaliblines(inccaliblines < (caliblines(end)+5));
    if nest
      newfolder = strcat(origpath,'\Calibration');
      newrawfldr = strcat(newfolder,'\RawData');
    else
      newfolder = strcat(origpath,'_calibration');
      newrawfldr = strcat(newfolder,'\RawData');
    end
    
    if ~isdir(newfolder)
      mkdir(newfolder);
    end
    
    if ~isdir(newrawfldr)
      mkdir(newrawfldr)
    end
    
    rawdatalines = unique(([caliblines;inccaliblines]));
    files2move = rawdatacontent(csvorder(rawdatalines));
    nlines = length(rawdatalines);
    newnums = cellstr(num2str([1:nlines]'));
    newnames = cell(nlines,1);
    newnames(:) = cellstr(repmat(rawhdrstr,nlines,1));
    newlinenames = cellfun(@(x,y) strcat(newrawfldr,'\',x,strtok(y),rawdatafmt),newnames,newnums,'uniformoutput',0);
    cellfun(@(x,y) copyfile(x,y),files2move,newlinenames)
    
    idcalib = find(cellfun(@(x) ~isempty(x),strfind(lower(line_types),'calib')));%ismember(lower(line_types),{'calib','solution'})));% 
    idnist = find(cellfun(@(x) ~isempty(x),strfind(lower(line_types),'nist')));
  elseif ~isempty(sollines)
    caliblines = sollines;
    nieghborlabels = unique(relhtxt([max(caliblines(1)-1,1),min(caliblines(end)+1,length(relhtxt))],1));
    inccaliblines = find(ismember(lower(relhtxt),'nist612'));
    inccaliblines = inccaliblines(inccaliblines < (caliblines(end)+5));
    if nest
      newfolder = strcat(origpath,'\Calibration');
      newrawfldr = strcat(newfolder,'\RawData');
    else
      newfolder = strcat(origpath,'_calibration');
      newrawfldr = strcat(newfolder,'\RawData');
    end
    
    if ~isdir(newfolder)
      mkdir(newfolder);
    end
    
    if ~isdir(newrawfldr)
      mkdir(newrawfldr)
    end
    
    rawdatalines = unique(([caliblines;inccaliblines]));
    files2move = rawdatacontent(csvorder(rawdatalines));
    nlines = length(rawdatalines);
    newnums = cellstr(num2str([1:nlines]'));
    newnames = cell(nlines,1);
    newnames(:) = cellstr(repmat(rawhdrstr,nlines,1));
    newlinenames = cellfun(@(x,y) strcat(newrawfldr,'\',x,strtok(y),rawdatafmt),newnames,newnums,'uniformoutput',0);
    cellfun(@(x,y) copyfile(x,y),files2move,newlinenames)
    
    idcalib = find(cellfun(@(x) ~isempty(x),strfind(line_types,'solution')));%ismember(lower(line_types),{'calib','solution'})));% 
    idnist = find(cellfun(@(x) ~isempty(x),strfind(lower(line_types),'nist')));
  else
    idcalib = [];
  end
  fldrpaths_relLA = cell(length(orderlines),1);
  
  lineoffset = 0;
  for k = orderlines(:)'
    if ismember(k,idcalib)
      if nest
        fldrpaths_relLA{k} = strcat('Calibration');
      else
        fldrpaths_relLA{k} = strcat(datasetfield,'_calibration');
      end
      %
      %       if ~isdir(newfolder)
      %         mkdir(newfolder);
      %       end
      %
      %       if ~isdir(newrawfldr)
      %         mkdir(newrawfldr)
      %       end
      %
      %       rawdatalines = find(linesind == k);
      %       files2move = rawdatacontent(csvorder(rawdatalines));
      %       nlines = length(rawdatalines);
      %       newnums = cellstr(num2str([1:nlines]' + lineoffset));
      %       newnames = cell(nlines,1);
      %       newnames(:) = cellstr(repmat(rawhdrstr,nlines,1));
      %       newlinenames = cellfun(@(x,y) strcat(newrawfldr,'\',x,strtok(y),rawdatafmt),newnames,newnums,'uniformoutput',0);
      %       cellfun(@(x,y) copyfile(x,y),files2move,newlinenames)
      %       lineoffset = lineoffset + nlines;
    else
      if nest
        newfolder = strcat(origpath,'\',line_types{k});
        newrawfldr = strcat(newfolder,'\RawData');
        fldrpaths_relLA{k} = strcat(line_types{k});
      else
        newfolder = strcat(origpath,'_',line_types{k});
        newrawfldr = strcat(newfolder,'\RawData');
        fldrpaths_relLA{k} = strcat(datasetfield,'_',line_types{k});
      end
      
      if ~isdir(newfolder)
        mkdir(newfolder);
      end
      
      if ~isdir(newrawfldr)
        mkdir(newrawfldr)
      end
      
      rawdatalines = find(linesind == k);
      files2move = rawdatacontent(csvorder(rawdatalines));
      nlines = length(rawdatalines);
      newnums = cellstr(num2str([1:nlines]'));
      newnames = cell(nlines,1);
      newnames(:) = cellstr(repmat(rawhdrstr,nlines,1));
      newlinenames = cellfun(@(x,y) strcat(newrawfldr,'\',x,strtok(y),rawdatafmt),newnames,newnums,'uniformoutput',0);
      cellfun(@(x,y) copyfile(x,y),files2move,newlinenames)
    end
  end
  
  if runread
    basefldr = unique(fldrpaths_relLA);
    usePropnames = ismember(PropertyNames,'comp') + ...
      ismember(PropertyNames,'dat') + ...
      ismember(PropertyNames,'parse') + ...
      ismember(PropertyNames,'line_offset');
    pnamesi = [1:2:(2*length(usePropnames));2:2:(2*length(usePropnames))];
    vinind = pnamesi(:,find(usePropnames == 0));
    vinind = vinind(:);
    for k = 1:length(basefldr)
      disp(sprintf('Writing to %s.........',basefldr{k}));
      if nest
        [d,hdrtxt] = readfiles('comp',strcat(origpath,'\'),...
          'dat',basefldr{k},'parse',0,varargin{vinind});
      else
        [d,hdrtxt] = readfiles('dat',basefldr{k},'parse',0,varargin{vinind});
      end
      
      if strmatch('plot',PropertyNames)
        close all;
        plot_heatmap(d,[],'thresh','auto',...
          'addlabels',hdrtxt,'save_heat',1);
      end
    end
    
  end
else
  disp('Only one label type present in the dataset. Data parsing is unnecessary');
end

end