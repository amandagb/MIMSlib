function xls2lzs(varargin)
%% Script: xls2lzs('file','filepath','Y',#,'X',#,'spotsize',#,'scanspeed',#,'linespacing',#)
% Description:  This script has two main functionalists:
%  
% 1) Select a file and the script will compute the run time for you based on the parameters specified in the file
% RUN IN COMMAND PROMPT: xls2lzs
% OUTPUT (TEXT PRINTED IN COMMAND WINDOW):
% Image Dimensions: 5980.00 um in y, 1000.00 um in x; 300 lines total
%  Spot size: 10 um, Scan speed: 10 um/sec, Line Spacing: 20 um
%  Time per line: 100.00 seconds, Time for total area:   8.33 hours
% 
% 2) User defined parameters- does not require an input file
% RUN IN COMMAND PROMPT: xls2lzs('file','none','Y',5980,'X',1000,'spotsize',5,'scanspeed',15,'linespacing',10)
% OUTPUT (TEXT PRINTED IN COMMAND WINDOW):
% Image Dimensions: 5980.00 um in y, 1000.00 um in x; 598 lines total
%  Spot size: 5 um, Scan speed: 15 um/sec, Line Spacing: 10 um
%  Time per line:  66.67 seconds, Time for total area:  11.07 hours
% 
% Example: xls2lzs('file','none','Y',5980,'X',1000,'spotsize',5,'scanspeed',15,'linespacing',10)
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   'filepath':  string indicating the excel document path which contains
%   the laser ablation method map. User may also input 'none', as the
%   PropertyValue if they wish to forego coversion of a file
%   'X', 'Y': integer indicating width (in x direction) and height (in y
%   direction) of LA area specified, in microns.  If filepath ~= 'none',
%   these will automatically be defined from the file specified. Otherwise,
%   DEFAULT will be Y = 6000, X = 1000 if use does not specify these
%   properties.
%   'scanspeed': integer indicating scan speed in microns per second
%   [DEFUALT = 10 um/sec]
%   'spotsize': integer indicating spot size in microns [DEFUALT = 10um]
%   'linespacing': integer indicating spacing between lines (center to
%   center spacing) in microns [DEFUALT = 2*spotsize]
%   'linetrantime': integer indicating the line transition time in seconds
%   [DEFUALT = 10]
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author            E-mail                      Version
%  10 Nov  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     2 - Added
%  run-time analysis

%% ---- Read data from excel files ----
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));
f = 'D:\My Documents\Dropbox\MADLab Research\Data\LA Setup Template\CODE EXAMPLE.xls';
if strmatch('file',PropertyNames)
  f = PropertyVal{strmatch('file',PropertyNames)};
else
  [f,fpath] = uigetfile('.xls','Select the directory containing the LA Method MAP you wish to convert to .lzs'); % Select directory with .Fin2 files
  f = strcat(fpath,f);
end

% while ~isdir(f)
%   f = uigetfile(current_dir,'Select the directory containing the LA Method MAP you wish to convert to .lzx'); % Select directory with .Fin2 files
% end

if ~strcmp('none',f)
  flzsnm = strcat(f(1:end-4),'TEST.lzs');
  flzs = fopen(flzsnm,'w+');
  [num,txt,raw] = xlsread(f);
  
  headers = raw(1,:);
  
  lzsheaders = {'ScanName','X1','Y1','Z1','X2','Y2','Z2','Method','SpotSize',...
    'SpaceBetweenSpots','SpaceBetweenLines','Energy','PulseRepRate','ScanRate',...
    'NumberOfShots','Defocus','Defocusvalue','HeliumFlowRate','GasBlank',...
    'PauseBetweenSamples','ShutterDelay','TriggerDelay', 'SampleRunTime',...
    'TotalSampleTime','NumberOfRuns','SegmentCount','SegmentNumber'};
  
  ordxlsh = {'Scan Name','X1','Y1','Z1','X2','Y2','Z2','Method',...
    'Aperture Size','Space Between Spots','Space Between Lines','Energy',...
    'Pulse Rep Rate','Scan Rate','Number of Shots','Defocus','Defocus Amount',...
    'He Flow Rate','Gas Blank','Pause Between Samples','Shutter Delay',...
    'Trigger Delay','Sample Run Time','Total Sample Time','Number of Runs'};
  
  x2lheadermap = cell(length(lzsheaders),3);
  for j = 1:length(headers)
    x2lheadermap(j,:) = {lzsheaders{j},ordxlsh{j},strmatch(ordxlsh{j},headers,'exact')};
  end
  j = 1 + j; x2lheadermap(j,1) = {lzsheaders{j}};
  j = 1 + j; x2lheadermap(j,1) = {lzsheaders{j}};
  
  fprintf(flzs,'<Sequence app="DigiLaz_III.exe" Version="1" SampleCount="300">\r\n');
  
  v = zeros(1,size(x2lheadermap,1));
  s = 1;
  i = 2;
  
  while s
    for h = 1:length(x2lheadermap)
      lzsfield = x2lheadermap{h,1};
      xlsind = x2lheadermap{h,3};
      if ~isempty(xlsind) % xlsind is empty for SegmentCount and SegmentNumber fields which both equal zero for all row, indicating that the values in v do not need to be adjusted
        if strmatch(lzsfield,'ScanName')
          if isnan(raw{i,xlsind})
            scanname = '';
          else
            scanname = raw{i,xlsind};
          end
          
        elseif strmatch(lzsfield,'Method')
          if strmatch(raw{i,xlsind},'Single Line');
            v(h) = 0;
          else
            v(h) = 0;
          end
          
        elseif ~isempty(strmatch(lzsfield,'SpaceBetweenSpots'))
          if isnan(raw{i,xlsind});
            xlsind = x2lheadermap{strmatch('SpotSize',x2lheadermap(:,1)),3};
            v(h) = raw{i,xlsind};
          else
            v(h) = raw{i,xlsind};
          end
          
        elseif ~isempty(strmatch(lzsfield,'SpaceBetweenLines'))
          if isnan(raw{i,xlsind});
            xlsind = x2lheadermap{strmatch('SpotSize',x2lheadermap(:,1)),3};
            v(h) = 2*raw{i,xlsind};
          else
            v(h) = raw{i,xlsind};
          end
          
        elseif strmatch(lzsfield,'NumberOfShots')
          if isnan(raw{i,xlsind});
            v(h) = 1;
          else
            v(h) = raw{i,xlsind};
          end
          
        elseif strmatch(lzsfield,'Defocus')
          if strmatch('N',raw{i,xlsind})
            v(h) = 0;
          else
            v(h) = 1;
          end
          
        else
          if isnan(raw{i,xlsind}) && ~isnan(v(h))
            s = 0;
          else
            v(h) = raw{i,xlsind};
          end
        end
      end
    end
    fprintf(flzs,'  <Row_%d ScanName="%s" X1="%6.2f" Y1="%6.2f" Z1="%6.2f" X2="%6.2f" Y2="%6.2f" Z2="%6.2f" Method="%d" SpotSize="%d" SpaceBetweenSpots="%d" SpaceBetweenLines="%d" Energy="%d" PulseRepRate="%d" ScanRate="%6.2f" NumberOfShots="%d" Defocus="%d" DefocusValue="%6.2f" HeliumFlowRate="%d" GasBlank="%d" PauseBetweenSamples="%d" ShutterDelay="%d" TriggerDelay="%d" SampleRunTime="%d" TotalSampleTime="%d" NumberOfRuns="%d" SegmentCount="%d" SegmentNumber="%d"/>\r\n',...
      i-1,scanname,v(2:end));
    i = i+1;
  end
  fclose(flzs);
else
  if strmatch('Y',PropertyNames)
    H = PropertyVal{strmatch('Y',PropertyNames)};
  else
    H = 6000;
  end
  
  if strmatch('X',PropertyNames)
    W = PropertyVal{strmatch('X',PropertyNames)};
  else
    W = 1000;
  end
  
  if strmatch('spotsize',PropertyNames)
    spotsize = PropertyVal{strmatch('spotsize',PropertyNames)};
  else
    spotsize = 10;
  end
  
  if strmatch('scanspeed',PropertyNames)
    scanspeed = PropertyVal{strmatch('scanspeed',PropertyNames)};
  else
    scanspeed = 10;
  end
  
  if strmatch('linespacing',PropertyNames)
    spacebtwnlns = PropertyVal{strmatch('linespacing',PropertyNames)};
  else
    spacebtwnlns = 2*spotsize;
  end
  nlines = floor(H/spacebtwnlns);
  R = 0;
end

%% Runtime Analysis
if ~strcmp('none',f)
  xlsind = find(strcmp('X1',headers) == 1);   X1 = cell2mat(raw(2:i-2,xlsind));
  xlsind = find(strcmp('X2',headers) == 1);   X2 = cell2mat(raw(2:i-2,xlsind));
  
  xlsind = find(strcmp('Y1',headers) == 1);   Y1 = cell2mat(raw(2:i-2,xlsind));
  xlsind = find(strcmp('Y2',headers) == 1);   Y2 = cell2mat(raw(2:i-2,xlsind));
  
  xlsind = find(strcmp('Aperture Size',headers) == 1);
  spotsize = cell2mat(raw(2,xlsind));     spacebtwnspots = spotsize;
  
  xlsind = find(strcmp('Scan Rate',headers) == 1);
  scanspeed = cell2mat(raw(2,xlsind));  nlines = length(X1);
  
  if Y1(1)-Y2(1) % indicates that X1(1) = X2(1), i.e., we are rastering in the horizontal (x) direction, for each line x is constant and y varies
    H = abs(Y1(1)-Y2(1));
    W = abs(X1(1)-X1(end));
    R = 1; % rastering in horizontal direction
    spacebtwnlns = abs(X1(1)-X1(2));
  elseif X1(1)-X2(1) % indicates that Y1(1) = Y2(1), i.e., we are rastering in the verticle (y) direction, for each line y is constant and x varies
    H = abs(Y1(1)-Y1(end)); % height in microns
    W = abs(X1(1)-X2(1)); % width in microns
    R = 0; % rastering in verticle direction
    spacebtwnlns = abs(Y1(1)-Y1(2));
  end
end

if strmatch('linetranstime',PropertyNames)
  linetranstime = PropertyVal{strmatch('linetranstime',PropertyNames)};
else
  linetranstime = 10; %[seconds]
end

if R
  tperline = H/scanspeed;
  ttotal = nlines*tperline + linetranstime;
else
  tperline = W/scanspeed;
  ttotal = nlines*tperline + linetranstime;
end

area_prop = sprintf('Image Dimensions: %6.2f um in y, %6.2f um in x; %6.2f lines total\n',H,W,nlines);
LA_params = sprintf('Spot size: %d um, Scan speed: %d um/sec, Line Spacing: %d um, Line Transition Time: %d sec\n',...
  spotsize,scanspeed,spacebtwnlns,linetranstime);
t_ana = sprintf('Time per line: %6.2f seconds, Time for total area: %6.2f hours',tperline,ttotal/3600);
disp(sprintf('%s %s %s',area_prop,LA_params,t_ana));
