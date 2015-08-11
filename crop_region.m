function [reginfo, fixed_reg ,moving_reg] = crop_region(photo,datasetnum,varargin)
%% Script: everything_out = EC720reg(datasetnum,varargin)
% [reginfo, fixed_reg ,moving_reg] = EC720reg([],1,'filtered',0,'element','Cu3273',...
%     'fix',1,'resize',1,'plot',1)
% center_imtrans_alpha(moving_image, [tx ty theta gx gy sx sy],...
%   'extrapval',#,'fixedstr',{'tx'},'fixedval','fixedvals',[#],...
%   'varstr',{'theta'},'varval',[5*pi/180])
% Description: This function wil take multiple images as inputs (were data
% inputs can be individual images or structures) and finds a common capture
% region which describes all the images
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% photo: either 2 or 3-D matrix corresponding to the photographic image. It
% could also be a string corresponding to the path & location of the image.
% datasetnum: number or string associated with the metallomic dataset
% varargin - 'PropertyName','PropertyValue'
%   'filtered': binary 0 or 1 to apply 'fix_sngl_pxl'
%   'element': string associated with the element that will be used
%   [DEFAULT = 'Cu3273']
%   'fix': binary indicating whether to fix photo (0) or MSI (1) [DEFAULT = 1]
%   'resize': binary indicating whether to resize images to match
%   dimensions of photo (0) or MSI (1) [DEFAULT = 1]
%   'plot': binary indicating whether to plot (1) or not (0) [DEFAULT = 1]
%   'exp#': integer indicating the experiment number
%   'fnum': integer indicating the experiment number
%   'mnum': integer indicating the experiment number
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  24 May 2012    Amanda Gaudreau   amanda.gaudreau@gmail.com     1

calibrate = 0;
if calibrate; load('D:\My Documents\MADLab Data\Laser Ablation\20120510_OESstandards\Analysis\calib_params.mat'); end;
% 20 um spot experiements
tags = {'20120320_OESCugrid20umps', '20120320_OESCugrid40umps',...
  '20120322_OESCugrid0,25ta', '20120320_OESCugrid100umps',...
  '20120322_OESCugrid1ta', '20120322_OESCugrid1ta3'};
nsets = length(tags);
el2ana = {'Si2516'};%{'Cu3247', 'Cu3273', 'Si2516'};
ind_el2ana = interp_data_elemrng(readfilesOES('comp','pc','dat',tags{1}),el2ana);
thresh = cell(1,length(el2ana));
data = cell(1,nsets);
medline = cell(1,nsets);
line = cell(1,nsets);
grid_region = zeros(nsets,2);
ncolumns = zeros(nsets,1);
for t = 1:nsets
  d = readfilesOES('comp','pc','dat',tags{t});
  data{t} = getfield(d,el2ana{1});
  ncolumns(t) = size(data{t},2);
  hline = sum(data{t});
  ord = round(size(data{t},1)*0.05);
  if rem(ord,2); ord = ord + 1; end;
  medline = medfilt1(hline,9);
  
  for j = 1:2
    dmedline = diff(medline);
    signmedline = sign(dmedline);
    %figure; plot(sum(data{t})); hold on; plot(medline{t},'r');
    [v,ind] = max(dmedline);
    climb = 1;
    strike = 0;
    while climb && strike < 3
      ind = ind + 1;
      if signmedline(ind) < 0
        strike = strike + 1;
      else
        strike = 0;
      end
    end
    if j == 1; hstart_ind = ind - strike; end
    if j == 2; hfinal_ind = length(medline) - (ind - strike) + 1; end
    medline = fliplr(medline);
  end
  hgrid_region(t,:) = [hstart_ind,hfinal_ind];
  width = hgrid_region(t,2) - hgrid_region(t,1);
  %plot_heatmap(d,el2ana,'scale',0,'plot','lin','thresh','auto');
  
end
width = hgrid_region(:,2) - hgrid_region(:,1);
nsidepnts = [hgrid_region(:,1) - 1,ncolumns - hgrid_region(:,2)];
perc = zeros(nsets,1);
perc_sides = zeros(nsets,2);
for t = 1:nsets
  perc_sides(t,:) = [nsidepnts(t,1)/width(t),nsidepnts(t,2)/width(t)];
  perc(t) = width(t)/ncolumns(t);
end
side_perc = floor(min([5,perc_sides(:)'*100]));
npnts = round(width*side_perc./100);

for t = 1:nsets
  d = readfilesOES('comp','pc','dat',tags{t});
  def_region = hgrid_region(t,1) - npnts(t):hgrid_region(t,2) + npnts(t);
%   data{t} = data{t}(:,def_region);
  plot_heatmap(d,el2ana,'scale',0,'plot','lin','thresh','auto','samples',def_region);
end
end