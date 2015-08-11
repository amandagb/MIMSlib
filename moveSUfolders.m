% function missingdata = renameivis(datestr,renamestr)
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

cdir = loaddirfun;
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
ivissucol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. spectral unmixed IVIS folder')));
fldrSUname = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Renamed spectral unmixed IVIS folder')));
ivisSUseq = XLtxt(:,ivissucol);
hSUseq = cellstr(strtok(ivisSUseq,'_'));

if exist('datestr') ~= 1
  datestr = '';
end

if isempty(datestr)
  alldates = unique(XLnum(:,datecol));
  datestr = num2str(unique(alldates));
elseif isnumeric(datestr)
  datestr = num2str(datestr);
end
datestr = cellstr(datestr);

if exist('renamestr') ~= 1
  renamestr = 'both';
end

[nrows,ncols] = size(XLall);
errstruct = struct;
for d = 1:length(datestr)
  datemat = str2num(datestr{d});
  ystr = datestr{d}(1:4);
  mstr = datestr{d}(5:6);
  dstr = datestr{d}(7:8);
  origivisfldrs = cellstr(ls(origivisdatadir));
  fldrstr = 'Spectral Unmixing Results for';
  keepfldrs = cellfun(@(x) ~isempty(x),strfind(origivisfldrs,fldrstr));
  origivisfldrs = origivisfldrs(keepfldrs);
  ymatch = cellfun(@(x) ismember(length(fldrstr)+1 + 1,x),strfind(origivisfldrs,ystr));
  mmatch = cellfun(@(x) ismember(length(fldrstr)+1 + 6,x),strfind(origivisfldrs,mstr));
  dmatch = cellfun(@(x) ismember(length(fldrstr)+1 + 9,x),strfind(origivisfldrs,dstr));
  fldrind = find(sum([ymatch,mmatch,dmatch],2) == 3);
  
  if ~isempty(fldrind)
    cd(origivisdatadir)
    cd(origivisfldrs{fldrind})
    SUchildren = cellstr(ls);
    SUchildren = SUchildren(3:end);
    seqfldrs = find(cellfun(@(x) ~isempty(x),strfind(SUchildren,'SEQ')));
    for f = seqfldrs(:)'
      cd(SUchildren{f})
      SUseqhdr = strtok(SUchildren{f},'_');
      SUdir = cd;
      seqChildren = cellstr(ls);
      seqinfofile = find(cellfun(@(x) ~isempty(x) && x == 1,strfind(seqChildren,'SequenceInfo.txt')));
      if ~isempty(seqinfofile)
        movefile(strcat(SUdir,'\SequenceInfo.txt'),strcat(SUdir,'\',SUseqhdr,'_SequenceInfo.txt'));
      end
      keepfldrs = cellfun(@(x) isempty(x),strfind(seqChildren,'.'));
      seqFolders = seqChildren(keepfldrs);
      hFolder = unique(strtok(seqFolders,'_'));
      if length(hFolder) > 1
        disp(sprintf('%s in %s has folders from %d sequences',...
          SUseqhdr,origivisfldrs{fldrind}(length(fldrstr)+1:end),length(hFolder)))
        seq2ivisfldr = cell(length(hFolder),2);
        for h = 1:length(hFolder)
          SUrow = find(cellfun(@(x) ~isempty(x),strfind(hSUseq,hFolder{h})));
          seq2ivisfldr{h,1} = seqChildren{find(cellfun(@(x) ~isempty(x),strfind(seqChildren,hFolder{h})))};
          if ~isempty(SUrow)
            seq2ivisfldr{h,2} = XLtxt(SUrow,fldrSUname);
          end
        end
        errstruct = setfield(errstruct,SUseqhdr,seq2ivisfldr);
        cd ..
      else
        SUrow = find(cellfun(@(x) ~isempty(x),strfind(hSUseq,hFolder{1})));
        if ~isempty(SUrow)
          ivisfldrname = XLtxt(SUrow,fldrSUname);
          ivisfldrpath = strcat(cdir.ivisfldr,ivisfldrname{1});
          movefile(strcat(SUdir,'\*'),ivisfldrpath);
          cd ..
          rmdir(SUdir)
        else
          cd ..
        end
      end
    end
  end
end
% end
