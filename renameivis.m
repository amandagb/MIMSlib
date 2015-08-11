function missingdata = renameivis(datestr,renamestr)
%% function: renameivis
% Example:
%   >> missingdata = renameivis;  -- runs through all dates in the excel
%           file and all possible data types (rgb, ivissu, ivisbs)
%   >> renameivis('20141118','rgb')
%
% Description: Uses information in the "IVISinfo" excel database to rename
% files and move them into a single folder.
%   The script REQUIRES the following folder structure:
%     YYYY-MM-DD IVIS and Gross Pathology --- folder containing the following folders:
%       YYYY-MM-DD IVIS data --- folder containing MADYYYYMMDDHHMMSS_SEQ folders
%       YYYY-MM-DD IVIS Gross Pathology --- folder containing #.TIF files
%     **** NOTE: all file names are case and white space sensitive (for now)
% The user must indicate a datestr and a target rename file type. The
% script uses the information within the IVISinfo.xls file to rename the
% files appropriately within the original folder structure. If the user
% indicates <<renamestr = 'rgb'>>, the corresponding tif file number must
% be indicated in the "Orig. TIF image number" column (currently col 13).
% If the user indicates <<renamestr = 'ivissu'>>, the sequence folder
% indicates in "Orig. spectral unmixed IVIS folder" column (col 17) must be
% indicated. In both cases, the files will be renamed according to the
% strings indicated in "RGB image name" and "Renamed spectral unmixed IVIS
% folder" respectively
%
% Required Functions: 
%     IVISinfo.xlsx ***** THIS FILE SHOULD RESIDE IN YOUR IVIS EXPERIMENT
%                         FOLDER. The script assumes this and uses the
%                         defined or indicated path as the file's path
%
% INPUTS ----------------------------------------------------------------
% The user should change the following varaibles to reflect their paths and
% desired analysis
% • datestr: string input with format 'YYYYMMDD'
% • renamestr: string indicating whether to rename tif images in "IVIS
%   Gross Pathology" folder ( = 'rgb') or IVIS sequence folders in the
%   "IVIS data" folder (='ivis' / 'ivissu' / 'ivisbs'). To do both, set
%   string to 'both'
% >>>>>> internal variables
% • deleteolddata: logical (1 or 0) indicating whether to delete old files
%   and directories after they have been moved
% • verbose:  logical (1 or 0) indicating whether error messages in
%   subfunctions should be broadcast to Command Window
% • IVISexpdir: this line may be uncommented and completed with the
%   appropriate path on users own PC
% • origivisdatadir: this line may be uncommented and completed with the
%   appropriate path on users own PC
% OUTPUTS ---------------------------------------------------------------
% • missingdata: structure with 6 fields (two per relevant data class --
% rgb, ivis spectral unmixing, ivis blue shift). Each field is a cell
% containing information about the data that are missing or are not
% defined. 
%   ----DefXL: indicates that these data are lacking definition in the
%   excel database file
%   ----Files: indicates that these files or folders do not exist in their
%   appropriate parent directories.
%  
%   COPIES indicated tif files to new directories under the IVISexpdir
%   **NOTE: these files can be moved rather than copied by changing the
%   "copyfile" command to "movefile"
%
%  Date           Author            E-mail                      Version
%  31 Jan  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      0
%   2 Feb  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      1
%     Changed so that if input is empty, all folders in the "ServerFiles"
%     directory are used
%  22 Feb  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      2
%     Adds spectral unmixed data to appropriate folders

deleteolddata = 1;
verbose = 1;
IVISexpdir = 'C:\Users\ADGB\Documents\MADLab Data\EB_Experiments';
if ~exist('IVISexpdir') || isempty(IVISexpdir)
  IVISexpdir = uigetdir(cd,'Select IVIS experiment parent directory'); % Select directory with .Fin2 files
end

origivisdatadir = 'C:\Users\ADGB\Documents\MADLab Data\EB_Experiments\ServerFiles';
if ~exist('origivisdatadir') || isempty(origivisdatadir)
  origivisdatadir = uigetdir(cd,'Select IVIS experiment archieve directory (may be the same as IVIS experiment directory)'); % Select directory with .Fin2 files
end

allrgbfldr = 'RGBmacroscopic';
allivisfldr = 'IVIS_EBFluor';
[XLnum,XLtxt,XLall] = xlsread(strcat(IVISexpdir,'\IVISinfo.xlsx'));
Xhdrs = XLall(1,:);

datecol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Sacrifice date')));
injurycol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Injury')));
expidcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Experiment ID')));
animalnumcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Animal ID number')));
anatomydircol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Anatomy depiction')));
injectioncol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Injection Type')));
origtifnumcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. TIF image number')));
polcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Polarizer')));
ivissucol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. spectral unmixed IVIS folder')));
ivisbscol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. blue shift IVIS folder')));
rgbimgcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'RGB')));

if exist('datestr') ~= 1
  datestr = '';
end

if isempty(datestr)
  alldates = unique(XLnum(:,datecol));
  datestr = num2str(unique(alldates));
elseif ~isstr(datestr)
  datestr = num2str(datestr);
end
datestr = cellstr(datestr);

if exist('renamestr') ~= 1
  renamestr = 'both';
end

[nrows,ncols] = size(XLall);

missingrgbdef = cell(nrows,ncols); % XLall(1:2,:);% 
missingrgbdef(1:2,:) = XLall(1:2,:);
missingrgbfiles = cell(nrows,ncols); % XLall(1:2,:);% 
missingrgbfiles(1:2,:) = XLall(1:2,:);
startrgbi = [3,3];
missingivissudef = cell(nrows,ncols); % XLall(1:2,:);% 
missingivissudef(1:2,:) = XLall(1:2,:);
missingivissufldr = cell(nrows,ncols); % XLall(1:2,:);% 
missingivissufldr(1:2,:) = XLall(1:2,:);
startivissui = [3,3];
missingivisbsdef = cell(nrows,ncols); % XLall(1:2,:);% 
missingivisbsdef(1:2,:) = XLall(1:2,:);
missingivisbsfldr = cell(nrows,ncols); % XLall(1:2,:);% 
missingivisbsfldr(1:2,:) = XLall(1:2,:);
startivisbsi = [3,3];

for d = 1:length(datestr)
  datemat = str2num(datestr{d});
  ystr = datestr{d}(1:4);
  mstr = datestr{d}(5:6);
  dstr = datestr{d}(7:8);
  origivisfldrs = cellstr(ls(origivisdatadir));
  keepfldrs = cellfun(@(x) length(x),strfind(origivisfldrs,'-')) == 2;
  origivisfldrs = origivisfldrs(keepfldrs);
  ymatch = cellfun(@(x) ismember(1,x),strfind(origivisfldrs,ystr));
  mmatch = cellfun(@(x) ismember(6,x),strfind(origivisfldrs,mstr));
  dmatch = cellfun(@(x) ismember(9,x),strfind(origivisfldrs,dstr));
  fldrind = find(sum([ymatch,mmatch,dmatch],2) == 3);
  
  if isempty(fldrind)
    renamefldrs = {''};
  else
    renamedir = strcat(origivisdatadir,'\',origivisfldrs{fldrind});
    renamefldrs = cellstr(ls(renamedir));
  end
  
  rgbfldrind = cellfun(@(x) ~isempty(x),strfind(renamefldrs,'Pathology'));
  if any(rgbfldrind)
    rgbfldr = strcat(renamedir,'\',renamefldrs{rgbfldrind});
  else
    rgbfldr = [];
    disp(sprintf('No %s-%s-%s IVIS Gross Pathology folder',ystr,mstr,dstr));
  end
  
  ivisfldrind = cellfun(@(x) ~isempty(x),strfind(renamefldrs,'data'));
  if any(ivisfldrind)
    ivisfldr = strcat(renamedir,'\',renamefldrs{ivisfldrind});
  else
    ivisfldr = [];
    disp(sprintf('No %s-%s-%s IVIS data folder',ystr,mstr,dstr));
  end
  
  switch renamestr
    case 'rgb'
      renamecell = {'rgb'};
    case 'ivissu'
      renamecell = {'ivissu'};
    case 'ivisbs'
      renamecell = {'ivisbs'};
    case 'ivis'
      renamecell = {'ivissu','ivisbs'};
    case 'both'
      renamecell = {'rgb','ivissu','ivisbs'};
    otherwise
      renamecell = {'rgb','ivissu'};
  end
  
  for f = 1:length(renamecell)
    renameftype = renamecell{f};
    
    xlsdates = XLnum(1:end,datecol);
    dateoffest = sum(cellfun(@(x) ~isempty(x),XLtxt(1:end,datecol)));
    daterows = find(xlsdates == datemat);
    
    switch renameftype
      %.....................................................................
      case 'rgb'
        tifnum = cellstr(num2str(XLnum(daterows,origtifnumcol)));
        tifnum = cellfun(@(x) strtok(x),tifnum,'uniformoutput',0);
        if any(rgbfldrind)
          newnames = XLtxt(daterows+dateoffest,rgbimgcol);
          [missingxlrgb,missingrgbf] = ...
            renamesubfun(renameftype,tifnum,newnames,Xhdrs,XLall(daterows+dateoffest,:),...
            IVISexpdir, allrgbfldr, renamedir,  renamefldrs{rgbfldrind},deleteolddata);
        else
          missingxlrgb = cellfun(@(x) ~isempty(x),strfind(tifnum,'NaN'));
          missingrgbf = true(length(daterows),1);
        end
        endi1 = startrgbi(1)+sum(missingxlrgb)-1;
        endi2 = startrgbi(2)+sum(missingrgbf)-1;
        missingrgbdef(startrgbi(1):endi1,:) = XLall(daterows(missingxlrgb)+dateoffest,:);
        missingrgbfiles(startrgbi(2):endi2,:) = XLall(daterows(missingrgbf)+dateoffest,:);
        startrgbi = [endi1+1,endi2+1];
        %.....................................................................
      case 'ivissu'
        ivisseq = XLtxt(daterows+dateoffest,ivissucol);
        if any(ivisfldrind)
          newnames = XLtxt(daterows+dateoffest,ivissucol+1);
          [missingxlivissu,missingivissuf] = ...
            renamesubfun(renameftype,ivisseq,newnames,Xhdrs,XLall(daterows+dateoffest,:),...
            IVISexpdir, allivisfldr, renamedir,  renamefldrs{ivisfldrind},deleteolddata);
        else
          missingxlivissu = cellfun(@(x) ~isempty(strmatch('',x)),ivisseq);
          missingivissuf = true(length(daterows),1);
        end
        endi1 = startivissui(1)+sum(missingxlivissu)-1;
        endi2 = startivissui(2)+sum(missingivissuf)-1;
        missingivissudef(startivissui(1):endi1,:) = XLall(daterows(missingxlivissu)+dateoffest,:);
        missingivissufldr(startivissui(2):endi2,:) = XLall(daterows(missingivissuf)+dateoffest,:);
        startivissui = [endi1+1,endi2+1];
        %.....................................................................
      case 'ivisbs'
        ivisseq = XLtxt(daterows+dateoffest,ivisbscol);
        if any(ivisfldrind)
          newnames = XLtxt(daterows+dateoffest,ivisbscol+1);
          [missingxlivisbs,missingivisbsf] = ...
            renamesubfun(renameftype,ivisseq,newnames,Xhdrs,XLall(daterows+dateoffest,:),...
            IVISexpdir, allivisfldr, renamedir,  renamefldrs{ivisfldrind},deleteolddata);
        else
          missingxlivisbs = cellfun(@(x) ~isempty(strmatch('',x)),ivisseq);
          missingivisbsf = true(length(daterows),1);
        end
        endi1 = startivisbsi(1)+sum(missingxlivisbs)-1;
        endi2 = startivisbsi(2)+sum(missingivisbsf)-1;
        missingivisbsdef(startivisbsi(1):endi1,:) = XLall(daterows(missingxlivisbs)+dateoffest,:);
        missingivisbsfldr(startivisbsi(2):endi2,:) = XLall(daterows(missingivisbsf)+dateoffest,:);
        startivisbsi = [endi1+1,endi2+1];
    end
  end
end
missingdata.rgbDefXL = missingrgbdef(1:startrgbi(1)-1,:);
missingdata.rgbFiles = missingrgbfiles(1:startrgbi(2)-1,:);
missingdata.IVISsuDefXL = missingivissudef(1:startivissui(1)-1,:);
missingdata.IVISsuFolders = missingivissufldr(1:startivissui(2)-1,:);
missingdata.IVISbsDefXL = missingivisbsdef(1:startivisbsi(1)-1,:);
missingdata.IVISbsFolders = missingivisbsfldr(1:startivisbsi(2)-1,:);


  function [missingXLdef,missingdata] = renamesubfun(indstr,fieldmatch,newnames,xlhead,xldata,...
      expdir,newfldr,origdir,origfldr,deleteold)
    
    subfdatecol = cellfun( @(x) ~isempty(x),strfind(xlhead,'Sacrifice date'));
    subfinjurycol = cellfun( @(x) ~isempty(x),strfind(xlhead,'Injury'));
    subfexpidcol = cellfun( @(x) ~isempty(x),strfind(xlhead,'Experiment ID'));
    subfanimalnumcol = cellfun( @(x) ~isempty(x),strfind(xlhead,'Animal ID number'));
    subfanatomydircol = cellfun( @(x) ~isempty(x),strfind(xlhead,'Anatomy depiction'));
    
    if ~exist('deleteold')
      deleteold = 0;
    end
    
    switch indstr
      case 'rgb'
        targetType = 'file';
        targetDes = 'RGB';
        targetEnd = '.tif';
        targetExtra = '';
      otherwise
        targetType = 'folder';
        targetDes = 'IVIS data';
        targetEnd = '\';
        targetExtra = '*';
    end
    
    if isempty(origfldr)
      disp(sprintf('No %s %ss available to rename',targetDes,targetType))
    else
      fullorigpath = strcat(origdir,'\',origfldr);
      newpath = strcat(expdir,'\',newfldr);
      
      if ~isdir(newpath)
        mkdir(expdir,newfldr);
      end
      missingXLdef = false(length(fieldmatch),1);
      missingdata = false(length(fieldmatch),1);
      
      for j = 1:length(fieldmatch)
        if ~isempty(fieldmatch{j}) && isempty(strmatch(fieldmatch{j},'NaN'))
          oldfullstr = sprintf('%s\\%s%s',fullorigpath,fieldmatch{j},targetEnd);
          newfullstr = sprintf('%s\\%s%s',newpath,newnames{j},targetEnd);
          if exist(oldfullstr) && ~exist(newfullstr)
            copyfile(sprintf('%s%s',oldfullstr,targetExtra),...
              newfullstr)
            if deleteold
              if isdir(oldfullstr)
                rmdir(oldfullstr,'s')
              elseif ~isdir(oldfullstr)
                delete(oldfullstr)
              end
            end
          elseif ~exist(oldfullstr) && exist(newfullstr)
            if verbose; disp(sprintf('A %s by the name of %s%s already exists',targetType,newnames{j},targetEnd)); end
          elseif ~exist(oldfullstr) && ~exist(newfullstr)
            if verbose; disp(sprintf('The %s %s (to be renamed %s) is missing from %s',targetType,fieldmatch{j},newnames{j},origfldr)); end
            missingdata(j) = 1;
          elseif exist(oldfullstr) && exist(newfullstr) && deleteold
            if verbose; disp(sprintf('A %s by the name of %s%s already exists, %s will be deleted',targetType,newnames{j},targetEnd,fieldmatch{j})); end
            if isdir(oldfullstr)
              rmdir(oldfullstr,'s')
            elseif ~isdir(oldfullstr)
              delete(oldfullstr)
            end
          elseif exist(oldfullstr) && exist(newfullstr) && ~deleteold
            if verbose; disp(sprintf('%s has already been moved and renamed to %s',fieldmatch{j},newnames{j})); end
          else
            xlinfo = xldata(j,:);
            if verbose; disp(sprintf('Unknwon error for %d_%s_%s%d%s',...
                xlinfo{subfdatecol}, xlinfo{subfinjurycol},...
                xlinfo{subfexpidcol},...
                xlinfo{subfanimalnumcol},xlinfo{subfanatomydircol})); end
          end
        else
          missingXLdef(j) = 1;
          xlinfo = xldata(j,:);
          if verbose; disp(sprintf('No %s %s is indicated for %d_%s_%s%d%s',...
              indstr,targetType,xlinfo{subfdatecol}, xlinfo{subfinjurycol},...
              xlinfo{subfexpidcol},...
              xlinfo{subfanimalnumcol},xlinfo{subfanatomydircol})); end
          %disp(sprintf('No original %s name (or number) is indicated (j = %d)',targetType,j))
        end
      end
    end
  end
end
