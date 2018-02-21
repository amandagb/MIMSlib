function [hdrfields,hdridx] = ICAPQhdrtxtfields

%  hdrfields: 
%  Date           Author            E-mail                      Version
%   4 Nov  2015   Amanda Balderrama amandagbalderrama@gmail.com     0

hdrfields = {'Line Name','Time Stamp','Name',...1-3
  'Version','File Version','Machine','Additional Gas Flow 1',...4-7
  'Additional Gas Flow 2','Additional Gas Flow 3','CCT Entry Lens',...8-10
  'Angular Deflection','Deflection Entry Lens','Extraction Lens 1 Polarity',...11-13
  'Extraction Lens 1 Negative','Extraction Lens 1 Positive','Spray Chamber Temperature',...14-16
  'Peristaltic Pump Speed','Cool Flow','Sampling Depth',...17-19
  'Plasma Power','Auxilliary Flow','Nebulizer Flow',...20-22
  'Torch Horizontal Position','Torch Vertical Position','Extraction Lens 2',...23-25
  'CCT Focus Lens','CCT Bias','CCT Exit Lens',...26-28
  'Focus Lens','D1 Lens','D2 Lens',...29-31
  'Quad Entry Lens','Pole Bias','CCT1 Flow',...32-34
  'CCT1 Shut-Off Valve','Virtual CCT Mass to Dac Factor','Virtual CCT Mass to Dac Offset',...35-37
  'Virtual CCT Mass parameter b','Virtual CCT Mass Maximum Dac Limit Set','RF Plasma Lit Readback',...38-40
  'RF FET Temperature Ok Readback','Plasma Power Readback','Pole Bias Readback',...41-43
  'Torch Horizontal Position Readback','Torch Vertical Position Readback','Sampling Depth Readback',...44-46
  'Extraction Lens 1 Negative Readback','Extraction Lens 1 Positive Readback','Extraction Lens 2 Readback',...47-49
  'Angular Deflection Readback','Quad Entry Lens Readback','Focus Lens Readback',...50-52
  'D2 Lens Readback','Deflection Exit Lens Readback','Deflection Entry Lens Readback',...53-55
  'D1 Lens Readback','CCT Exit Lens Readback','CCT Bias Readback',...56-58
  'CCT Entry Lens Readback','CCT Focus Lens Readback','Analyzer Vacuum Ok Readback',...59-61
  'Interface Pressure Readback','Analyzer Pressure Readback','Detector Voltage (Counting) Readback',...62-64
  'Detector Voltage (Analog) Readback','Plasma Cooling Water Flow Readback','Interface Temperature Readback',...65-67
  'Exhaust Flow Readback','Quad Board Temperature Readback','Inlet Fan Speed Readback',...68-70
  'Outlet Fan Speed Readback','Peltier Temperature Hot Side Readback','Spray Chamber Temperature Readback',...71-73
  'Supply Voltage 500 V Readback','Supply Voltage 1 kV Readback','Nebulizer Supply Pressure Readback',...74-76
  'Nebulizer Flow Readback','Cool Flow Readback','Auxilliary Flow Readback',...77-79
  'Gas Flow 2 Readback','CCT1 Flow Readback','Additional Gas Flow 1 Readback',...80-82
  'Threshold','Element Names','dwell time','xcal factor'};%83-86
hdridx.LineName = find(cellfun(@(x) ~isempty(x) && x == 1,strfind(hdrfields,'Line Name')));
hdridx.TimeStamp = find(cellfun(@(x) ~isempty(x) && x == 1,strfind(hdrfields,'Time Stamp')));
hdridx.MfctrName = find(cellfun(@(x) ~isempty(x) && x == 1,strfind(hdrfields,'Name')));
hdridx.InstStr =  find(cellfun(@(x) ~isempty(x) && x == 1,strfind(hdrfields,'Machine')));
end