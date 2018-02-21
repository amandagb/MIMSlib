function [Sampdata,SampSum,EXdata,ExSum] = extractXLinfo(ExpID,varargin)
%% Script: [Sampdata,SampSum,EXdata,ExSum] = extractXLinfo(ExpID,varargin)
% Description: Reads the excel file containing inforation related to the
% sring "ExpID". As of the file version, there are 4 ExpID possibilities:
% 'MIMS', 'ALON', 'LacSac' and 'MMM'. The Samp--- output structures contain
% information on the "SampleSummary" sheet and the EX--- output structures
% contain information on the MIMS shet. The excel files should contain
% information important for spatial scaling, sample screening, experimental
% information, etc. 
% 
% Example: 
% 
% INPUTS ----------------------------------------------------------------
% ExpID:  string indicating the experiment ID you'd like to read the excel
%         database of. Current options include:
%   • MIMS
%   • ALON
%   • LacSac
%   • MMM
% varargin - 'PropertyName','PropertyValue'
%   • 'MIMS': EXdata sturcture outputted by extractXLinfo
%   • 'Samp': Sampdata sturcture outputted by extractXLinfo
%   • 'seqfldr':  string containing the folder of the sequence you'd like
%   the row information of
% 
% NOTE: these are found the "MADLab Data" parent directory with a YYYYMM_ExpID 
%   formatting assumed. The excel file should be INSIDE the experiment's
%   parent directory and the name should also contain the ExpID string and
%   the word "Log" in the file name.
%
% OUTPUTS ---------------------------------------------------------------
% Sampdata: structure with the fields {num,txt,raw} containing the raw
%   outputs from the xlsread function for the "SampleSummary" sheet
% SampSum: structure with fields containing the information for the
%   SampleSummary sheet. Field values are integers corresponding to the row
%   or column 
%     Fields names are as follows:
% EXdata: structure with the fields {num,txt,raw} containing the raw
%   outputs from the xlsread function for the "MIMS" sheet
% ExSum: structure with fields containing the information for the MIMS
%   sheet. Field values are integers corresponding to the row or column
% 
%
%  Date           Author            E-mail                      Version
%  3  June 2016   Amanda Balderrama amandagbalderrama@gmail.com     0

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

cdir = loaddirfun;
mdatals = cellstr(ls(cdir.mdatapath));
mdatafldrs = mdatals(find(cellfun(@(x) isempty(strfind(x,'.')),mdatals)));
Expfldr = mdatafldrs( cellfun(@(x) ~isempty(strfind(x,ExpID)),mdatafldrs));
Expls = cellstr(ls(strcat(cdir.mdatapath,Expfldr{1},'\')));
Expfile = Expls(cellfun(@(x) ~isempty(strfind(lower(x),'log')),Expls));
fname = strcat(cdir.mdatapath,Expfldr{1},'\',Expfile{1});

if strmatch('MIMS',PropertyNames)
  EXdata = PropertyVal{strmatch('MIMS',PropertyNames)};
  ndataMIMS = EXdata.num;
  txtMIMS = EXdata.txt;
  alldatMIMS = EXdata.raw;
else
  [ndatMIMS, txtMIMS, alldatMIMS] = xlsread(fname,'MIMS');
  EXdata.num = ndatMIMS;
  EXdata.txt = txtMIMS;
  EXdata.raw = alldatMIMS;
end

if strmatch('Samp',PropertyNames)
  Sampdata = PropertyVal{strmatch('Samp',PropertyNames)};
  ndatSamp = Sampdata.num;
  txtSamp = Sampdata.txt;
  dataSamp = Sampdata.raw;
else
  [ndatSamp, txtSamp, dataSamp] = xlsread(fname,'SampleSummary');
  Sampdata.num = ndatSamp;
  Sampdata.txt = txtSamp;
  Sampdata.raw = dataSamp;
end

[hdrrow,rotcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Section Orientation')));
[~,injcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Injury')));
[~,sampnum] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Sample')));
[~,seqfldrcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Sequence Folder')));
[~,tposcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Type Position')));
[~,spotcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Spot size')));
[~,vscancol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Scan Speed')));
[~,dl2lcol] = find(cellfun(@(x) ~isempty(x),strfind(txtSamp,'Line to line')));

SampSum.headerRow = hdrrow;
SampSum.rotation = rotcol;
SampSum.injection = injcol;
SampSum.sampnum = sampnum;
SampSum.seqfldr = seqfldrcol;
SampSum.typenum = tposcol;
SampSum.spotsize = spotcol;
SampSum.vscan = vscancol;
SampSum.dl2l = dl2lcol;

[~,newlables] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'New Sample Labels')));
[labrow,dsetcol] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Dataset Name')));
[~,isocol] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Isotope List')));
[~,calibcol] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Standard Concentration Values (ppm)')));
[~,useablecol] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Usable?')));

ExSum.headrRow = labrow;
ExSum.seqfldr = dsetcol;
ExSum.isolist = isocol;
ExSum.uselogical = useablecol;
ExSum.calibvals = calibcol;
[~,ExSum.origlabels] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Original Sample Labels')));
[~,ExSum.newlables] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'New Sample Labels')));
[~,ExSum.numelements] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'# elements')));
[~,ExSum.isodwellT] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Dwell Times')));
[~,ExSum.totaqhours] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'Total Hours')));
% [~,ExSum.] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'')));
% [~,ExSum.] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'')));
% [~,ExSum.] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'')));
% [~,ExSum.] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'')));
% [~,ExSum.] = find(cellfun(@(x) ~isempty(x),strfind(txtMIMS,'')));

samplabels = txtMIMS(labrow+2:end,newlables);
isostrs = txtMIMS(labrow+2:end,isocol);
alliso = unique(strsplit( cell2mat(cellfun(@(x) strrep(strcat(x,','),' ',''),isostrs,'uniformoutput',0)') ,','));
isocell = cellfun(@(x) strtrim(strsplit(x,',')),isostrs,'uniformoutput',0);
sampcommas = strfind(samplabels,',');
typesCell = cellfun(@(x) strsplit(x,','),samplabels,'uniformoutput',0);

EXdata.seqfldrs = txtMIMS(labrow+2:end,dsetcol);
EXdata.isostrs = isostrs;
EXdata.uniqueiso = alliso;
EXdata.typelabels = samplabels;
EXdata.typeCell = typesCell;

if strmatch('seqfldr',PropertyNames)
  dataset = PropertyVal{strmatch('seq',PropertyNames)};
  dlog = cellfun(@(x) ~isempty(strfind(x,dataset)),ASdata.txt);
  [dr,dc] = find(dlog);
  AnSum.datasetRow = dr;
  
  
end


end