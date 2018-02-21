function [dout,varargout] = rotMIMS(data,rotTag)
%% Script:  data_struct = rotMIMS(data,rotTag);
% Description: Rotates the MIMS data structure according to the roatation
% tag (defined visually in "MIMS Experiment Summary and Notes.pdf")
%
% Example:
%
% Required Function:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% rotTag: rotation "label" indicating how data need to be rotated s.t.
%       the anatomical left is on the the viewer's left and anatomy are
%       oriented dorsal on top to ventral on bottom (that is NoEaSoWe =>
%       Dorsal, Right, Ventral, Left)
%       -2    Original data aquired VRDL, flip columns up/down (or flip LR
%             then rotate twice)
%       -1    Original data aquired RVLD, rotate 90 deg clockwise
%        0    Original data aquired DRVL
%        1    Original data aquired LDRV, rotate 90 deg ccw
%        2    Original data aquired VLDR, rotate 180 deg
% varargin - 'PropertyName','PropertyValue'
%   • 'labelstr': string indicating the corresponding label map to isol
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author            E-mail                      Version
%  31 July 2015   Amanda Balderrama amanda.gaudreau@gmail.com     0

if isstruct(data) 
  dout = rmfield(data,skip_fields);
  dout.Time = data.Time;
  switch rotTag
    case -2
      dout = structfun(@(x) rot90(fliplr(x),rotTag),dout);
    otherwise
      dout = structfun(@(x) rot90(x,rotTag),dout);
  end
  dout.line_datevec = data.line_datevec;
  dout.dataset = data.dataset;
elseif iscell(data)
  switch rotTag
    case -2
      dout = cellfun(@(x) rot90(fliplr(x),rotTag),data);
    otherwise
      dout = cellfun(@(x) rot90(x,rotTag),data);
  end
else
  nch = size(data,3);
  dout = zeros([size(rot90(data(:,:,1),rotTag)),nch]);
  for c = 1:nch
    switch rotTag
      case -2
        dout(:,:,c) = rot90(fliplr(data(:,:,c)),rotTag);
      otherwise
        dout(:,:,c) = rot90(data(:,:,c),rotTag);
    end
  end
end
