function lzsout = readlzsfiles(lzsfile,varargin)
%% Script:
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   'filepath':  string indicating the excel document path which contains
%   the laser ablation method map
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author              E-mail                        Version
%  16 Jan  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
% lzsfile = '2015-01-05 MIMS 37 lzs.lzs';

compdir = loaddirfun;

if isstruct(lzsfile)
  LAdir = strcat(compdir.lafldr,lzsfile.dataset,'\RawData\');
  lzsfile = strcat(lzsfile.dataset,'.lzs');
elseif strmatch('\',lzsfile)
  LAdir = lzsfile;
  lzsfile = '';
else
  LAdir = strcat(compdir.lafldr,'lzsfiles\');
end

%% ---- Read data from excel files ----
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

A = importdata(strcat(LAdir,lzsfile));
Hdr = A{1};
A = A(2:end);
nArows = size(A,1);

lzsheaders = {'Row','ScanName=','X1=','Y1=','Z1=','X2=','Y2=','Z2=','Method=','SpotSize=',...
  'SpaceBetweenSpots=','SpaceBetweenLines=','Energy=','PulseRepRate=','ScanRate=',...
  'NumberOfShots=','Defocus=','DefocusValue=','HeliumFlowRate=','GasBlank=',...
  'PauseBetweenSamples=','ShutterDelay=','TriggerDelay=', 'SampleRunTime=',...
  'TotalSampleTime=','NumberOfRuns=','SegmentCount=','SegmentNumber='};
nfields = length(lzsheaders);
% linefmt = strcat('<Row_%d ScanName="%s" X1="%1.2f" Y1="%1.2f" Z1="%1.2f"',...
%   ' X2="%1.2f" Y2="%1.2f" Z2="%1.2f" Method="%d" SpotSize="%d" ',...
%   ' SpaceBetweenSpots="%d" SpaceBetweenLines="%d" Energy="%d" ',...
%   ' PulseRepRate="%d" ScanRate="%1.2f" NumberOfShots="%d" Defocus="%d" ',...
%   ' DefocusValue="%1.2f" HeliumFlowRate="%d" GasBlank="%d" ',...
%   ' PauseBetweenSamples="%d" ShutterDelay="%d" TriggerDelay="%d" ',...
%   ' SampleRunTime="%d" TotalSampleTime="%d" NumberOfRuns="%d" ',...
%   ' SegmentCount="%d" SegmentNumber="%d"/>');
hlinefmt = '<Sequence app="%s" Version="%s" SampleCount="%s">';
starti = strfind(lower(Hdr),'samplecount')+length('samplecount')+2;
endi = strfind(Hdr(starti:end),'"');
nsamples = str2num(Hdr(starti:(starti + endi - 2)));

strind = nan(nArows,nfields+1);
for f = 1:nfields
  indout = cellfun(@(x) strfind(x,lzsheaders{f}) + length(lzsheaders{f}) + 1 ,A,'uniformoutput',0);
  indoutmtx = cell2mat(indout);
  nonemptycell = find(cell2mat(cellfun(@(x) ~isempty(x),indout,'uniformoutput',0)));
  strind(nonemptycell,f) = indoutmtx;
end
indout = cellfun(@(x) length(x),A,'uniformoutput',0);
indoutmtx = cell2mat(indout);
nonemptycell = find(cell2mat(cellfun(@(x) ~isempty(x),indout,'uniformoutput',0)));
strind(nonemptycell,end) = indoutmtx;
fnumbers = nan(nArows,nfields);
for r = 1:nArows
  Aline = A{r};
  for f = 1:nfields
    %Sf = cellfun(@(x) strfind(x(strind(f):strind(f+1)),'"'),A,'uniformoutput',0);
    startindf = strind(r,f);
    if ~isnan(startindf)
      if f > 1
        qind = strfind(Aline(startindf:end),'"');
        endindf = startindf + qind(1) - 2;
      else
        endindf = strind(r,f+1) - length(lzsheaders{f+1})-3;
      end
      Aftxt = Aline(startindf:endindf);
      Aval = str2num(Aftxt);
      if ~isempty(str2num(Aftxt))
        fnumbers(r,f) = Aval;
      end
    end
  end
end
keeprows = find(sum(isnan(fnumbers),2) < nfields);
keepflds = find(sum(isnan(fnumbers(keeprows,:)))==0);
fnumbers = fnumbers(keeprows,keepflds);
lzsheaders = lzsheaders(keepflds);

x1 = unique(fnumbers(:,strmatch('x1',lower(lzsheaders))));
y1 = unique(fnumbers(:,strmatch('y1',lower(lzsheaders))));
z1 = unique(fnumbers(:,strmatch('z1',lower(lzsheaders))));

x2 = unique(fnumbers(:,strmatch('x2',lower(lzsheaders))));
y2 = unique(fnumbers(:,strmatch('y2',lower(lzsheaders))));
z2 = unique(fnumbers(:,strmatch('z2',lower(lzsheaders))));

xyz = zeros(2,3);
if length(x1) > 1
  xyz(1,1) = x1(2) - x1(1);
  xyz(2,1) = x2(2) - x2(1);
else
  xyz(1,1) = x1;
  xyz(2,1) = x2;
end

if length(y1) > 1
  xyz(1,2) = y1(2) - y1(1);
  xyz(2,2) = y2(2) - y2(1);
else
  xyz(1,2) = y1;
  xyz(2,2) = y2;
end

if length(z1) > 1
  xyz(1,3) = z1(2) - z1(1);
  xyz(2,3) = z2(2) - z2(1);
else
  xyz(1,3) = z1;
  xyz(2,3) = z2;
end

spotsize = unique(fnumbers(:,strmatch('spotsize',lower(lzsheaders))));

spacebtwnspots = unique(fnumbers(:,strmatch('spacebetweenspots',lower(lzsheaders))));
spacebtwnlines = unique(fnumbers(:,strmatch('spacebetweenlines',lower(lzsheaders))));
energy = unique(fnumbers(:,strmatch('energy',lower(lzsheaders))));
scanspeed = unique(fnumbers(:,strmatch('scanrate',lower(lzsheaders))));

lzsout.headers = lzsheaders;
lzsout.data = fnumbers;
lzsout.xyz = xyz;
lzsout.spotsize = spotsize;
lzsout.spacebtwnspots = spacebtwnspots;
lzsout.spacebtwnlines = spacebtwnlines;
lzsout.energy = energy;
lzsout.scanspeed = scanspeed;
lzsout.filename = lzsfile;
end


