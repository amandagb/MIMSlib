function varargout = getphysicalparams(data,varargin)
%% Script:  paramStruct = getphysicalparams(d.dataset)
% Description:  Given the dataset tag, uses the information in
% CBM Data Log.xls to determine what the physical parameters of the
% experiment.
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure output by readfiles function
% varargin - 'PropertyName','PropertyValue'
%       CURRENTLY NONE
% OUTPUTS ---------------------------------------------------------------
% paramStruct: structure containing fields with relevant spatial
% information for the experiment
%   
%
%  Date           Author              E-mail                      Version
%  6  June 2012   Amanda Gaudreau     amanda.gaudreau@gmail.com       1
%  20 Jan  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%   9 Jan  2017   Amanda Balderrama   amandagbalderrama@gmail.com     3
%     Uses MIMSsummary.m file to search for parameters
%  13 Jan  2017   Amanda Balderrama   amandagbalderrama@gmail.com     3.1
%  24 Feb  2017   Amanda Balderrama   amandagbalderrama@gmail.com     3.2
%   Small adjustment to error prompt

cdir = loaddirfun;
mimsInfo = MIMSsummary;
pspotPos = 1;
pl2lPos = 2;
pscanPos = 3;
ppowerPos = 4;
pflowPos = 5;
psampPos = 6;

if isstruct(data)
  tag = data.dataset;
  T = data.Time;
elseif isstr(data)
  tag = data;
else
  tag = 'notag';
end
dpos = find(cellfun(@(x) ~isempty(strfind(x,tag)),mimsInfo.dsetInfo(:,1)));

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('outstruct',PropertyNames)
  outstruct = PropertyVal{strmatch('outstruct',PropertyNames)};
else
  outstruct = 1;
end

rotval = 0;
errflag = 0;
if isempty(dpos) && exist(strcat(cdir.mlabpath,'CBM Data Log v3.xls'),'file')
  [N,T,R] = xlsread(strcat(cdir.mlabpath,'CBM Data Log v3.xls'),'DATA LOG 2012');
  data_tags = T(5:end,2);
  header_offset = 4;
  tfcell = strcmp(data_tags,tag);
  tfvec = tfcell;
  %tfvec = cellfun(@(x) ~isempty(x),tfcell);
  dataind = find(tfvec == 1);
  if isempty(dataind)
    disp(sprintf('The dataset %s does not exist in the CBM Data Log excel file. Please double check the dataset''s name.',tag));
    parameters = [50,50,50,100,1.5,1];
    errflag = 1;
  else
    parameters = N(dataind,6:11); % (1) spot size, (2) l2l distance, (3) scan speed, (4) power, (5) gas flow rate, (6) sampling rate
  end
elseif ~isempty(dpos)
  parameters = nan(1,6);
  parameters(pspotPos) = mimsInfo.dsetInfo{dpos(1),mimsInfo.dspotCOL};
  parameters(pl2lPos) = mimsInfo.dsetInfo{dpos(1),mimsInfo.dl2lCOL};
  parameters(pscanPos) = mimsInfo.dsetInfo{dpos(1),mimsInfo.vscanCOL};
  parameters(psampPos) = nan;
  rotval = mimsInfo.rotvec(dpos(1));
else
  errflag = 1;
  parameters = [50,50,50,100,1.5,nan];
  disp({[sprintf('Acquisition parameters for %s are not defined',tag)];...
    ['in the "CBM Data Log v3.xls" or the "MIMSsummary.m" file,'];...
    ['default parameters are provided:'];...
    [sprintf('spot = %d um, scan = %d um/s, l2l = %d um',parameters([pspotPos,pscanPos,pl2lPos]))]});
end

if exist('T')
  if isempty(parameters(psampPos)) || isnan(parameters(psampPos))
    T(T == 0) = nan;
    Tdiff = diff(T')';
    samptime = abs(nanmean(Tdiff(:)));
    if samptime > 1e-3
      parameters(psampPos) = samptime;
    else
      parameters(psampPos) = nanmean(nanmean(diff(T)));
    end
  end
end

if nargout == 1;
  if outstruct
    pstruct.spotsize = parameters(1);
    pstruct.l2ldist = parameters(2);
    pstruct.scanspeed = parameters(3);
    pstruct.energy = parameters(4);
    pstruct.gasrate = parameters(5);
    pstruct.tsamp = parameters(6);
    pstruct.rotval = rotval;
    pstruct.errorflag = errflag;
    varargout{1} = pstruct;
  else
    varargout{1} = parameters;
  end
else
  i = 1;
  varargout{i} = parameters(1); i = i + 1; % (1) spot size
  varargout{i} = parameters(2); i = i + 1; % (2) l2l distance
  varargout{i} = parameters(3); i = i + 1; % (3) scan speed
  varargout{i} = parameters(4); i = i + 1; % (4) power
  varargout{i} = parameters(5); i = i + 1; % (5) gas flow rate
  varargout{i} = parameters(6); i = i + 1; % (6) sampling rate
  varargout{i} = errflag;
end
end