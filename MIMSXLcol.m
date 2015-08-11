function colstruct = MIMSXLcol(hdrtxt)
%% function: colstruct = ivisXLcol
%
% Description: Gives the column index in the headerinfo.xlsx file
%       'Line Name', 'Time Stamp', 'Name', 'Version', 'File Version',
%       'Machine', 'Additional Gas Flow 1', 'Additional Gas Flow 2',
%       'Additional Gas Flow 3', 'CCT Entry Lens', 'Angular Deflection',
%       'Deflection Entry Lens', 'Extraction Lens 1 Polarity', 'Extraction
%       Lens 1 Negative', 'Extraction Lens 1 Positive', 'Spray Chamber
%       Temperature', 'Peristaltic Pump Speed', 'Cool Flow', 'Sampling
%       Depth', 'Plasma Power', 'Auxilliary Flow', 'Nebulizer Flow', 'Torch
%       Horizontal Position', 'Torch Vertical Position', 'Extraction Lens
%       2', 'CCT Focus Lens', 'CCT Bias', 'CCT Exit Lens', 'Focus Lens',
%       'D1 Lens', 'D2 Lens', 'Quad Entry Lens', 'Pole Bias', 'CCT1 Flow',
%       'CCT1 Shut-Off Valve', 'Virtual CCT Mass to Dac Factor', 'Virtual
%       CCT Mass to Dac Offset', 'Virtual CCT Mass parameter b', 'Virtual
%       CCT Mass Maximum Dac Limit Set', 'RF Plasma Lit Readback', 'RF FET
%       Temperature Ok Readback', 'Plasma Power Readback', 'Pole Bias
%       Readback', 'Torch Horizontal Position Readback', 'Torch Vertical
%       Position Readback', 'Sampling Depth Readback', 'Extraction Lens 1
%       Negative Readback', 'Extraction Lens 1 Positive Readback',
%       'Extraction Lens 2 Readback', 'Angular Deflection Readback', 'Quad
%       Entry Lens Readback', 'Focus Lens Readback', 'D2 Lens Readback',
%       'Deflection Exit Lens Readback''Deflection Entry Lens Readback',
%       'CCT Exit Lens Readback', 'CCT Bias Readback', 'CCT Entry Lens
%       Readback', 'CCT Focus Lens Readback', 'Analyzer Vacuum Ok
%       Readback', 'Interface Pressure Readback', 'Analyzer Pressure
%       Readback', 'Detector Voltage (Counting) Readback', 'Detector
%       Voltage (Analog) Readback', 'Plasma Cooling Water Flow Readback',
%       'Interface Temperature Readback', 'Exhaust Flow Readback', 'Quad
%       Board Temperature Readback', 'Inlet Fan Speed Readback', 'Outlet
%       Fan Speed Readback', 'Peltier Temperature Hot Side Readback',
%       'Spray Chamber Temperature Readback', 'Supply Voltage 500 V
%       Readback', 'Supply Voltage 1 kV Readback', 'Nebulizer Supply
%       Pressure Readback', 'Nebulizer Flow Readback', 'Cool Flow
%       Readback', 'Auxilliary Flow Readback', 'Gas Flow 2 Readback', 'CCT1
%       Flow Readback', 'Additional Gas Flow 1 Readback', 'Threshold',
%       'Element Names', 'dwell time', 'xcal factor'
% Example:
% Required Functions: loaddirfun, IVISinfo.xlsx
%
% INPUTS ----------------------------------------------------------------
% 
% OUTPUTS ---------------------------------------------------------------
% colstruct
%
%  Date           Author            E-mail                      Version
%   1 Apr  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      0

compdir = loaddirfun;

if ~exist('hdrtxt')
  [~,~,hdrtxt] = xlsread('headerinfo.xlsx');
  Xhdrs = hdrtxt(1,:);
end

colstruct.datecol = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Time Stamp')));
colstruct.linelabel = find( cellfun( @(x) ~isempty(x),strfind(Xhdrs,'Line Name')));
end