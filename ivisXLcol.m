function colstruct = ivisXLcol(Xhdrs)
%% function: colstruct = ivisXLcol
%
% Description: Gives the column index in the IVISinfo.xlsx file
%   {'Sacrifice date','Injury','Experiment ID','Animal ID number','Anatomy
%   depiction','Injection Type','Injected Quantity (mL)','BUCS
%   Pre-score','BUCS Post1-score','BUCS Post2-score','BUCS
%   Post3-score','BUCS 3hr-score','Orig. TIF image number','Polarizer','RGB
%   image name','Polarization Type','Orig. spectral unmixed IVIS
%   folder','Renamed spectral unmixed IVIS folder','Orig. blue shift IVIS
%   folder','Renamed blue shift IVIS folder'}
% Example:
% Required Functions: loaddirfun, IVISinfo.xlsx
%
% INPUTS ----------------------------------------------------------------
% 
% OUTPUTS ---------------------------------------------------------------
% colstruct
%
%  Date           Author            E-mail                      Version
%  25 Feb  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      0

compdir = loaddirfun;

if ~exist('Xhdrs')
  [XLnum,XLtxt,XLall] = xlsread(strcat(compdir.EBexppath,'IVISinfo.xlsx'));
  Xhdrs = XLall(1,:);
end

colstruct.datecol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Sacrifice date')));
colstruct.injurycol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Injury')));
colstruct.include = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Include')));
colstruct.expidcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Experiment ID')));
colstruct.animalnumcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Animal ID number')));
colstruct.anatomydircol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Anatomy depiction')));
colstruct.injectioncol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Injection Type')));
colstruct.origtifnumcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. TIF image number')));
colstruct.polcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Polarizer')));
colstruct.bucscols = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'BUCS')));
colstruct.origivissuSEQ = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. spectral unmixed IVIS folder')));
colstruct.origivisbsSEQ = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Orig. blue shift IVIS folder')));
colstruct.rgbimgcol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'RGB')));
end