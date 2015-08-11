function [flds,varargout] = skip_fields(varargin)
%% Script: flds = skip_fields
% Description:  Catalogs the fields which should be skipped. In all files
% that cycle through an analysis of all element fields of a data structure,
% this function will be required. The fields are as follows and will be
% updated as new non-element fields are added to the data structure:  
%         {'Time','line_datevec','dataset'}
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   {1}:  structure generated from readfiles.m script which
%   contains fields of all elements collected from OES or ICP-MS and other
%   informative fields such as Time, line_datevec and dataset
% OUTPUTS ---------------------------------------------------------------
% flds: cell of strings with fields that should be skipped
% varargout
%   {1-total number of flds to skip}: if a data structure is also entered
%   as an input, the user should also specify outputs which will carry the
%   data from the fields that should be skipped
%
%  Date           Author            E-mail                      Version
%  3  July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

flds = {'Time','line_datevec','dataset'};

if ~isempty(varargin)
  dstruct = varargin{1};
  
  PropertyNames = varargin(2:2:length(varargin));
  PropertyVal = varargin(3:2:length(varargin));
  
  if strmatch('data',PropertyNames)
    nbins = PropertyVal{strmatch('data',PropertyNames)};
  else
    nbins = 1e3;
  end
  
  [elem_ana_range,fldnm] = interp_data_elemrng(data,flds);
  i = 1;
  for f = elem_ana_range
    varargout{i} = getfield(dstruct,fldnm{f});
    i = i + 1;
  end
end
end
