%% Script: reformatnames
%
% Description: Finds the BUCS score information for a given dataset
% identified by the header "animalquery". The script reads through the
% excel file "EBexpinfo.xlsx". The first two rows of the file are headers:
%   {'Sacrifice date','Animal ID number','Comments','Anatomy depiction',
%   'EB injection','EB quantity (ml/kg)','Injury','Pad',
%   'BUCS Pre-score (open, mesh, beam)','BUCS Post1-score (open, mesh, beam)',
%   'BUCS Post2-score (open, mesh, beam)','BUCS Post3-score (open, mesh, beam)',
%   'BUCS 3hr-score (open, mesh, beam)','RGB image name','Polarizer',
%   'Polerization Type','IVIS Fluorescent image sequence'}
% Example:
%     >> animalstr = '20140819_nopad5D';
%     >> [bucs,bstr,banimals] = extractBUCS(animalstr);
%     >> [bucs,bstr,banimals] = extractBUCS('L_');
% Required Functions: loaddirfun, IVISinfo.xlsx
%
% INPUTS ----------------------------------------------------------------
% varargin
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author            E-mail                      Version
%  31 Jan  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      0

IVISexpdir = 'C:\Users\ADGB\Documents\MADLab Data\EB_Experiments';
allrgbfldr = 'RGBmacroscopic';
allivisfldr = 'IVIS_EBFluor';
[XLnum,XLtxt,XLall] = xlsread(strcat(IVISexpdir,'\IVISinfo.xlsx'));
ndatacol = size(XLall,2);
ndatarows = size(XLtxt,1);
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
bucscol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'BUCS')));
rgbimgcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'RGB')));

xlsdates = XLnum(1:end,datecol);
dateoffest = sum(cellfun(@(x) ~isempty(x),XLtxt(1:end,datecol)));
xlsannum = XLnum(1:end,animalnumcol);
xlsanatomy = XLtxt(1:end,anatomydircol);
xlsinjury = XLtxt(1:end,injurycol);

injurymatch = {'pad','nopad','blast','b','sham','s','sb','si','nc'};
injurylength = cellfun(@(x) length(x),injurymatch);
anatomydir = {'D_','V_','L_','R_','S_'};
polstr = 'pol';
psstr = 'PS';

cd(strcat(IVISexpdir,'\',allrgbfldr))
rgbfiles = cellstr(ls);
keepfiles = cellfun(@(x) ~isempty(strfind(lower(x),'.tif')),rgbfiles);
rgbfiles = rgbfiles(keepfiles);
for f = 1:length(rgbfiles)
  fname = rgbfiles{f};
  fdate = str2num(fname(1:8));
  maskpos = strfind(fname,psstr);
  if isempty(maskpos)
    fnamebase = fname(1:end-4);
    maskstr = '';
  else
    fnamebase = fname(1:(maskpos-2));
    maskstr = fname((maskpos-1):(end-4));
  end
  fdirpos = cellfun(@(x) strfind(fnamebase,x),anatomydir,'UniformOutput',0);
  fdiruse = cellfun(@(x) ~isempty(x),fdirpos);
  fanatomydir = anatomydir{fdiruse}(1);
  %injurypos = cellfun(@(x) strfind(fnamebase,x),injurymatch,'UniformOutput',0);
  
  numendi = fdirpos{fdiruse}-1;
  numstarti = numendi;
  while ~isempty(str2num(fname(numstarti:numendi)))
    numstarti = numstarti - 1;
  end
  numstarti = numstarti + 1;
  animalid = str2num(fname(numstarti:numendi));
  daterows = find(xlsdates == fdate);
  annumrows = find(xlsannum == animalid);
  relrows = daterows(ismember(daterows,annumrows));
  filerow = relrows(ismember(xlsanatomy(relrows+dateoffest),fanatomydir));
  xlsanimalrow = XLall(filerow+dateoffest,:);
  xlsnamefmt = xlsanimalrow{rgbimgcol};
  
  if isempty(strfind(fname,xlsnamefmt))
    newname = strcat(xlsnamefmt,maskstr,'.tif');
    movefile(fname,newname)
    disp(newname);
  end
end

cd(strcat(IVISexpdir,'\',allivisfldr))
ivisfldrs = cellstr(ls);
keepfiles = cellfun(@(x) ~isempty(strfind(lower(x),'.tif')),rgbfiles);
rgbfiles = rgbfiles(keepfiles);
for f = 1:length(rgbfiles)
  fname = rgbfiles{f};
  fdate = str2num(fname(1:8));
  maskpos = strfind(fname,psstr);
  if isempty(maskpos)
    fnamebase = fname(1:end-4);
    maskstr = '';
  else
    fnamebase = fname(1:(maskpos-2));
    maskstr = fname((maskpos-1):(end-4));
  end
  uscorei = strfind(fname,'_');
  fdirpos = cellfun(@(x) strfind(fnamebase,x),anatomydir,'UniformOutput',0);
  fdiruse = cellfun(@(x) ~isempty(x),fdirpos);
  fanatomydir = anatomydir{fdiruse}(1);
  
  %injurypos = cellfun(@(x) strfind(fnamebase,x),injurymatch,'UniformOutput',0);
  
  numendi = fdirpos{fdiruse}-1;
  numstarti = numendi;
  while ~isempty(str2num(fname(numstarti:numendi)))
    numstarti = numstarti - 1;
  end
  numstarti = numstarti + 1;
  animalid = str2num(fname(numstarti:numendi));
  daterows = find(xlsdates == fdate);
  annumrows = find(xlsannum == animalid);
  relrows = daterows(ismember(daterows,annumrows));
  filerow = relrows(ismember(xlsanatomy(relrows+dateoffest),fanatomydir));
  xlsanimalrow = XLall(filerow+dateoffest,:);
  xlsnamefmt = xlsanimalrow{rgbimgcol};
  
  if isempty(strfind(fname,xlsnamefmt))
    newname = strcat(xlsnamefmt,maskstr,'.tif');
    movefile(fname,newname)
    disp(newname);
  end
end

