function [elem_ana_range,fldnm] = interp_data_elemrng(data,v)
%% Script:  [elem_ana_range,fldnm] = interp_data_elemrng(data,v)
% Description:  Used to read data from the data structure created in
% readfiles.m.
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% data: a structure contining a field with the name of each
% element for which data was collected. Each field contains a <# lines> x
% <# time samples> double matrix. (output from readfiles.m)
% v: 1 x <# elements> matrix which corresponds to the field number in the data structure. User may also specify
%   This variable can also be a cell of strings or a single string where the string specifies
%   the isotope the user would like to analyze.
% OUTPUTS ---------------------------------------------------------------
% elem_ana_range:
%
%  Date           Author            E-mail                      Version
%  19 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  10 Mar  2015   AmandaGBalderrama amanda.gaudreau@gmail.com     1.1 
if isstruct(data)
  fldnm = fieldnames(rmfield(data,skip_fields));
  elem_ana_range = [];
  if ~isnumeric(v)
    if iscell(v)
      for i = 1:length(v)
        if isempty(strmatch(lower(v{i}),lower(fldnm)))
          disp(['The element you have indicated, (', v{i},'), does not exist in the dataset']);
        else
          elem_ana_range = [elem_ana_range,strmatch(lower(v{i}),lower(fldnm))'];
        end
      end
    else
      elem_ana_range = [elem_ana_range,strmatch(lower(v),lower(fldnm))];
      if isempty(elem_ana_range)
        elem_ana_range = 1;
      end
    end
  elseif isempty(v)
    %COMMENTED OUT 15 SEPT 2011 for id_MIMS_regions
    %if isempty(strmatch('time',lower(fldnm)));
    %  elem_ana_range = 1;
    %else
      elem_ana_range = 1:length(fldnm);
    %end
  else
      elem_ana_range = v;
  end
else
  if isempty(v)
    fldnm{1} = 'din'; % input data
  else
    fldnm{1} = v;
  end
  elem_ana_range = 1;
end

elem_ana_range = elem_ana_range(:)';
end
