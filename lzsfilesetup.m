function xls2lzs(varargin)
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
%  Date           Author            E-mail                      Version
%  9  Nov  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

%% ---- Read data from excel files ----
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));
f = 'D:\My Documents\Dropbox\MADLab Research\Data\LA Setup Template\CODE EXAMPLE.xls';
if strmatch('file',PropertyNames)
  f = PropertyVal{strmatch('file',PropertyNames)};
else
  [f,fpath] = uigetfile('.xls','Select the directory containing the LA Method MAP you wish to convert to .lzx'); % Select directory with .Fin2 files
  f = strcat(fpath,f);
end

% while ~isdir(f)
%   f = uigetfile(current_dir,'Select the directory containing the LA Method MAP you wish to convert to .lzx'); % Select directory with .Fin2 files
% end

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
%% Runtime Analysis
xlsind = find(strcmp('X1',headers) == 1);   X1 = cell2mat(raw(2:i-2,xlsind));
xlsind = find(strcmp('X2',headers) == 1);   X2 = cell2mat(raw(2:i-2,xlsind));

xlsind = find(strcmp('Y1',headers) == 1);   Y1 = cell2mat(raw(2:i-2,xlsind));
xlsind = find(strcmp('Y2',headers) == 1);   Y2 = cell2mat(raw(2:i-2,xlsind));

xlsind = find(strcmp('Aperture Size',headers) == 1);   
spotsz = cell2mat(raw(2:i-2,xlsind));     spacebtwnspots = spotsz;    
spacebtwnlns = 2*spotsz;

xlsind = find(strcmp('Scan Rate',headers) == 1);   
scanspd = cell2mat(raw(2:i-2,xlsind));     

xlsind = find(strcmp('Pause Between Samples',headers) == 1);   
tbtwnpxls = cell2mat(raw(2:i-2,xlsind));     

if Y1(1)-Y2(1)
  
elseif X1(1)-X2(1)
  
end
