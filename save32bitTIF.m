function save32bitTIF(data_struct,data_str,varargin)
%% function:  save32bitTIF(data_struct,data_str)
% Description:  
% Example:  
% INPUTS ----------------------------------------------------------------
% data_struct:  a structure contining a field with the name of each
%   element for which data was collected. Each field contains a <# lines> x
%   <# time samples> double matrix.
% data_str: string modifier that is concatonated at the end of element
%   string name
% varargin - 'PropertyName','PropertyValue'
%   • 'saveto': string indicating the folder name that the images should be
%       saved to [DEFAULT = 'Images']
%   • 'dataheader': logical indicating whether to include the dataheader
%     string (from data_struct.dataset) to the beginning of the tif file
%     name (1) or not (0) [DEFAULT = 0]
% 
% OUTPUTS ---------------------------------------------------------------
% 
%  Date           Author            E-mail                      Version
%  12 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     0
%  13 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     0.1
%     cd-ed back to original directory
%  17 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     1
%     Allowed variable number of inputs for increased personalization
%   6 Mar  2015   Amanda Balderrama amanda.gaudreau@gmail.com     1.1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('saveto',PropertyNames)
  fldrname = PropertyVal{strmatch('saveto',PropertyNames)};
else
  fldrname = 'Images';
end

if strmatch('dataheader',PropertyNames)
  dataheader = PropertyVal{strmatch('dataheader',PropertyNames)};
else
  dataheader = 1;
end

if strmatch('basedir',PropertyNames)
  cdir.lafldr = PropertyVal{strmatch('basedir',PropertyNames)};
  if isempty(strmatch(cdir.lafldr(end),'\'))
      cdir.lafldr = strcat(cdir.lafldr,'\');
  end
else
  cdir = loaddirfun;
end

cd(cdir.lafldr);
cd(data_struct.dataset);
del = rmfield(data_struct,skip_fields);
delnms = fieldnames(del);
delcell = struct2cell(del);
[M,N] = size(delcell{1});
if ~isdir(strcat(cd,'\',fldrname))
  mkdir(strcat(cd,'\',fldrname));
end
cd(fldrname)
for i = 1:length(delnms)
  if dataheader
    tiffname = sprintf('%s %s%s.tif',data_struct.dataset,delnms{i},data_str);
  else
    tiffname = sprintf('%s%s.tif',delnms{i},data_str);
  end
  
  if exist(tiffname) ~= 2
    t = Tiff(tiffname,'w');
    t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
    t.setTag('Compression',Tiff.Compression.None);
    t.setTag('BitsPerSample',32);
    t.setTag('SamplesPerPixel',1);
    t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
    t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
    t.setTag('ImageLength',M);
    t.setTag('ImageWidth',N);
    %t.setTag('TileLength',32);
    %t.setTag('TileWidth',32);
    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    t.write(uint32(delcell{i}));
    t.close();
  end
end
cd ..