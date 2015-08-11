function [d,varargout] = cps2ppm(dcps,dppm,varargin)
%% Script: d = reshape_OESstructure(dstand,varargin)
% Description: Takes the OES standard text file with structure
% LABELline#,MM/DD/YYYY,
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% dcps:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed (+ other fields such as Time, dataset,
%       line_datevec)
% dppm: data structure with fields of elements for which calibration
%       parameters are available. Each dppm.ELEMENT will itself have two
%       fields: line_coeff and rsq
%         • line_coeff is the output from polyfit which takes as in input the
%         mean signal for each of the signals (in cps) and the absolute
%         concentration values (in ppm) --> line_coeff = polyfit(ppm', stnd_mean(:,f),1)
%         • rsq is the 
%         ** See MADLABmtng_20120501 lines 255 - 263 for details
% varargin - 'PropertyName','PropertyValue'
%   'plot': binary indicating whether to generate debugging plots (1) or
%   not (0) [DEFAULT = 0]
%   'save_data': binary indicating whether to save restructured data (1) or
%   not (0) [DEFAULT = 0]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  18 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

fldnm_cps = fieldnames(dcps);
fldnm_ppm = fieldnames(dppm);

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 1;
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 0;
end

d = [];
for f = 1:length(fldnm_cps)
  if isempty(strmatch(fldnm_cps{f},skip_fields)) && f <= length(fldnm_cps) 
    if strmatch(fldnm_cps{f},fldnm_ppm)
      data_cps = getfield(dcps,fldnm_cps{f});
      linparams = getfield(dppm,fldnm_cps{f});
      fitparams = getfield(linparams,'lin_coeff');
      a1 = fitparams(1);
      a0 = fitparams(2);
      p1 = 1/a1;
      p0 = -a0/a1;
      data_ppm = polyval([p1,p0],data_cps);
      d = setfield(d,fldnm_cps{f},data_ppm);
    else
      disp(sprintf('%s in cps data structure does not have ppm curve parameters',fldnm_cps{f}))
    end
  else
    data = getfield(dcps,fldnm_cps{f});
    d = setfield(d,fldnm_cps{f},data);
  end
end

if svdat
  home_dir = strcat('D:\My Documents\MADLab Data\Laser Ablation\',d.dataset,'\');
  allfiles = ls;
  files = cellstr(allfiles);
  fldnm = fieldnames(d);
  for i = 1:length(fldnm)
    if isempty(strmatch(fldnm{i},{'dataset'}))
      csvname = strcat(home_dir,fldnm{i},'_ppm.csv');
      csvwrite(csvname,getfield(d,fldnm{i}));
    end
  end
end

if p
  plot_heatmap(d,[],'scale',0,'plot','lin','thresh','auto','label','c');
end

end
