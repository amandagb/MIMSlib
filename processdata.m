function [d,varargout] = processdata(din,varargin)
%% Script: [d,varargout] = processdata(din,varargin)
% Description: Takes the OES standard text file with structure
% LABELline#,MM/DD/YYYY,
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'linespersample':
%   'temp_el': string indicate the template element to be used. The string
%   must uniquely identify the element and must have the first letter
%   capitalized. [DEFAULT = 'C_']
%   'plot': binary indicating whether to generate debugging plots (1) or
%   not (0) [DEFAULT = 0]
%   'save_heat': binary indicating whether to save heatmaps (1) or not (0)
%   [DEFAULT = 0]
%   'save_data': binary indicating whether to save restructured data (1) or
%   not (0) [DEFAULT = 0]
%   'eval_samp_num': evaluates specific sample numbers [DEFAULT = evaluate
%   all sample numbers]
%   'temp_lines': determine which lines to use for the template
%   'keep_lines': vector containing lines that should be kept of final
%   image [DEFAULT = 'all']
%   'nalign': indicates the number of points that are allowed to vary
%   between the rises of subsequent lines [DEFAULT = 20]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  14 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

d = readfiles('comp','pc','dat',data_tag,'readcsv',0,'save_data',0);
plot_heatmap(d,'Si','plot','lin','thresh','auto','scale',0);
ds = reshape_OEScorr(d,'temp_el','Si','lines',2,...'keep_lines',1:75,...
  'plot',0,'save_data',0,'save_heat',0,'top_perc',45,'bottom_perc',25);

