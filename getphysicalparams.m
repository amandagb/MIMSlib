function varargout = getphysicalparams(data,varargin)
%% Script:  [dspot, dl2l, vscan, E, vgas, ta] = getphysicalparams(d.dataset)
% Description:  Given the dataset tag, uses the information in
% CBM Data Log.xls to determine what the physical parameters of the
% experiment.
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% tag:  string indicating the dataset name (usually of format
% YYYYMMDD_INSTRtissuetype)
% varargin - 'PropertyName','PropertyValue'
%       CURRENTLY NONE
% OUTPUTS ---------------------------------------------------------------
% varargout: either one single variable which will be a 1x6 vector with
% parameter in the following order: (1) spot size, (2) l2l distance, (3)
% scan speed, (4) power, (5) gas flow rate, (6) acquisition time; ALSO can
% be a 6 variable output with each variable being one of the six parameter
% in the same order above.
%
%  Date           Author              E-mail                      Version
%  6  June 2012   Amanda Gaudreau     amanda.gaudreau@gmail.com       1
%  20 Jan  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2

cdir = loaddirfun;

if isstruct(data)
  tag = data.dataset;
  T = data.Time;
elseif isstr(data)
  tag = data;
end

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('outstruct',PropertyNames)
  outstruct = PropertyVal{strmatch('outstruct',PropertyNames)};
else
  outstruct = 0;
end

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

if exist('T')
  if isempty(parameters(6)) || isnan(parameters(6))
    Tdiff = diff(data.Time')';
    meancols = mean(Tdiff(:,1:end-1),1);
    meancols = meancols(isnan(meancols) == 0);
    meanrows = mean(Tdiff(:,1:end-1),2);
    meanrows = meanrows(isnan(meanrows) == 0);
    parameters(6) = mean([meancols(:);meanrows(:)]);
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