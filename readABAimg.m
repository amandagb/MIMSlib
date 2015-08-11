function [RGB_img] = readABAimg(varargin);
%% Script:  [RGB_img] = readABAimg(varargin);
% Description:  
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   'zoom':  numeric specifying zoom (corresponds to resolution) of img.
%   Lower numbers correspond to lower resolutions (ie lowest rez. for 0)
%   [DEFAULT = -1, highest res. img]
%   'imgstr':
%   'imgid':
%   'position':
% OUTPUTS ---------------------------------------------------------------
% 
%
%  Date           Author            E-mail                      Version
%  28 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  02 Aug  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1.1
%     New version not created. XML link updated to match new syntax. See

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

%% ISH and Expression data
ImageSeriesID = 69873899;  % the ImageSeriesID is associated with Slc30a3-Sagittal-05-1346 --> Zn transporter 3
xDocstr = strcat('http://api.brain-map.org/api/v2/section_image_download/',num2str(ImageSeriesID));%,'?downsample=[#]&quality=[#]');
  % OLD XML SYNTAX: strcat('http://www.brain-map.org/aba/api/imageseries/',num2str(ImageSeriesID),'.xml');
  % NEW XML SYNTAX: 'http://api.brain-map.org/api/v2/section_image_download/[SectionImage.id]?downsample=[#]&quality=[#]&expression=[true|false]'
xDoc = xml_read(xDocstr);

% Only need one of these image indicators, this should become an input!!!
imgstr = 'Slc30a3_57';
imgid = 69834530; % Corresponds to Slc30a3_57
position = 1625;

numimgs = length(xDoc.images.image);
imgstruct = xDoc.images.image;
imgind = 0;
if ~isempty(position)
  while ~imgind
    for i = 1:numimgs
      imgdata = imgstruct(i,1);
      if position == imgdata.position.CONTENT
        imgind = i;
      else clear imgdata
      end
    end
  end
elseif ~isempty(imgid)
  while ~imgind
    for i = 1:numimgs
      imgdata = imgstruct(i,1);
      if imgid == imgdata.imageseriesid.CONTENT
        imgind = i;
      else clear imgdata
      end
    end
  end
elseif ~isempty(imgstr)
  while ~imgind
    for i = 1:numimgs
      imgdata = imgstruct(i,1);
      if strcmpi(imgstr,imgdata.imagedisplayname)
        imgind = i;
      else clear imgdata
      end
    end
  end
end

path.ISH =strcat('&path=',imgstruct(imgind,1).downloadImagePath);%/external/aibssan/production1/Slc30a3_05-1346_835/zoomify/primary/0202023009/Slc30a3_57_0202023009_D.aff';
path.EXP =strcat('&path=',imgstruct(imgind,1).downloadExpressionPath);%/external/aibssan/production1/Slc30a3_05-1346_835/zoomify/primary/0202023009/Slc30a3_57_0202023009_D.aff';
if strmatch('z',PropertyNames)
  zoom = PropertyVal{strmatch('z',PropertyNames)};
else
  zoom = -1;
end
zoom = strcat('zoom=',num2str(zoom));
mime = '';  %mime = '&mime=2';
top = '';   %top = '&top=0';
left = '';  %left = '&left=0';
width = ''; %width = '&width=0';
height = '';  %height = '&height=0'
api_root = 'http://www.brain-map.org/aba/api/image?';
api_url.ISH = strcat(api_root,zoom,path.ISH,mime,top,left,width,height);
api_url.EXP = strcat(api_root,zoom,path.EXP,mime,top,left,width,height);

% Read in and save images from the given API URL
RGB_img.ISH = imread(api_url.ISH);
RGB_img.EXP = imread(api_url.EXP);

%Determine the correct pixel scaling given the zoom level
[r,c,d] = size(RGB_img.ISH);
topzoom_umperpxl = 1.07;  %microns per pixel of maximum size image;
scale = imgstruct(imgind,1).height.CONTENT/r;
umperpxl = topzoom_umperpxl*scale;
RGB_img.params.umperpxl = umperpxl;
RGB_img.params.imageseriesdisplay = xDoc.imageseriesdisplayname;
RGB_img.params.imageseriesid = xDoc.imageseriesid.CONTENT;


%% Reads in the closes corresponding ABA Reference Atlas image based on the
% position of the indicated image
xARstr = 'http://mouse.brain-map.org/atlas/ARA/Sagittal/browser.xml';  %AR = Atlas Ref
xAR = xml_read(xARstr);
ARid = imgstruct(imgind,1).referenceatlasindex.CONTENT;
nAR = length(xAR.atlasimage);
ARpos = [0:200:(nAR-1)*200];  % reference images depict slices of mouse brain every 200 um
d = abs(imgstruct(imgind,1).position.CONTENT-ARpos);
[v,ind] = min(d);
ARimgind = ind;
ARimgstruct = xAR.atlasimage;

path.AR =strcat('&path=',ARimgstruct(ARimgind,1).referenceatlasfilename);%/external/aibssan/production1/Slc30a3_05-1346_835/zoomify/primary/0202023009/Slc30a3_57_0202023009_D.aff';
api_url.AR = strcat(api_root,zoom,path.AR,mime,top,left,width,height);
RGB_img.AR = imread(api_url.AR);
end