function [filestruct,queryXL] = getfilelist(varargin)
%% function: filestruct = getfilelist(varargin)
%
% filestruct = getfilelist('anatomy','V','injury','pad','polarizer',2);
%
% Description: Finds information for a given dataset
% identified by the header "animalquery". The script reads through the
% excel file "IVISinfo.xlsx". The first two rows of the file are headers:
%   {'Sacrifice date','Injury','Experiment ID','Animal ID number','Anatomy
%   depiction','Injection Type','Injected Quantity (mL)','BUCS
%   Pre-score','BUCS Post1-score','BUCS Post2-score','BUCS
%   Post3-score','BUCS 3hr-score','Orig. TIF image number','Polarizer','RGB
%   image name','Polarization Type','Orig. spectral unmixed IVIS
%   folder','Renamed spectral unmixed IVIS folder','Orig. blue shift IVIS
%   folder','Renamed blue shift IVIS folder'}
% Required Functions: loaddirfun, IVISinfo.xlsx
%
% INPUTS ----------------------------------------------------------------
% varargin - strings with additional modifiers for the rgb image followed
% by 'PropertyName','PropertyValue' pairs
%     • 'date': NUMERIC data that is either a single number or a vector 
%     (1 x ndates) of numbers indicating which matching dates user would 
%     like to query
%         Ex: 'date',[20140819,20141212]
%     • 'injury': STRING indicating type of injury ==> ERRORS****
%         Ex: 'injury',{'pad','nopad','s','si','sb','b'}
%     • 'experiment': STRING indicating experiement identifier (case
%     sensitive)
%         Ex: 'experiment','IVIS'
%     • 'animal': NUMERIC (either scalar or vector) indicating animal
%     numers to query
%         Ex: 'animal',[21,22,23]
%     • 'anatomy': STRING or cell of strings indicating anatomical
%     directions to extract
%         Ex: 'anatomy',{'D','V','L','R','S'}
%     • 'bucs':  ==> ERRORS****
%     • 'polarizer': NUMERIC indicating polarizer type used
%         {'0: none; 1: adequate; 2: better; 3: best;'}
%         EX: 'polarizer',2
%           
% OUTPUTS ---------------------------------------------------------------
% filestruct: structure with two fields 
%     • tifFiles: cell with file names in the 'RGB image name' column of the
%       IVISinfo.xlsx file which match query input
%     • ivisSUfldr: cell with folder names in the 'Renamed spectral unmixed
%       IVIS folder' column of the IVISinfo.xlsx file which match query input
%
%  Date           Author            E-mail                      Version
%  22 Feb  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      1

% PropertyNames = {'anatomy','polarizer','injury'};
% PropertyVal = {'V',2,'pad'};
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

umxedind = strmatch('unmixed',PropertyNames);
if any(umxedind)
  PropertyNames{umxedind} = 'umxed';
end

cdir = loaddirfun;
IVISexpdir = cdir.EBexppath;
allrgbfldr = 'RGBmacroscopic';
allivisfldr = 'IVIS_EBFluor';
[XLnum,XLtxt,XLall] = xlsread(strcat(IVISexpdir,'IVISinfo.xlsx'));
ndatacol = size(XLall,2);
ndatarows = size(XLtxt,1);
Xhdrs = cellstr(XLall(1,:));
[M,N] = size(XLall);
numcol = false(1,N);
strcol = true(1,N);
ncol = find(sum(isnan(XLnum)) ~= size(XLnum,1));
numcol(ncol) = true;
strcol(ncol) = false;

filetifcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'RGB image name')));
fldrivissucol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Renamed spectral unmixed IVIS folder')));

hdroffset = M - size(XLnum,1);
valmatch = zeros(M-hdroffset,1);
npropmatch = 0;
for i = 1:length(PropertyNames)
  hdrcol = cellfun(@(x) ~isempty(x),strfind(lower(Xhdrs),PropertyNames{i}));
  qfield = XLall((hdroffset+1):end,hdrcol);
  if numcol(hdrcol)
    struse = cellfun(@(x) ismember(x,PropertyVal{i}),qfield);
    npropmatch = npropmatch + 1;
    valmatch = struse + valmatch;
  elseif strcol(hdrcol)
    strval = cellfun(@(x) strfind(x,PropertyVal{i}) == 1,qfield,'uniformoutput',0);
    struse = cellfun(@(x) ~isempty(x) && x == 1,strval);
    npropmatch = npropmatch + 1;
    valmatch = struse + valmatch;
  end
end
relrows = find(valmatch == npropmatch) + hdroffset;
filestruct.tifFiles = XLall(relrows,filetifcol);
filestruct.ivisSUfldr = XLall(relrows,fldrivissucol);
queryXL = XLall([1,2,relrows(:)'],:);
% filestruct.ivisBSfldr = XLall(relrows,fldrivissucol);

end