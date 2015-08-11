function [Calstruct,dcellOut,varargout] = getCalibrationFit(d_div,line_types,laONmasks,varargin)
%% function: getCalibrationFit
% Description:
% Example:
%
%
% INPUTS ----------------------------------------------------------------
% datastruct:    MIMS structure OR a string indiciating the dataset the
% user would like to use.
% linelabels:
% varargin - 'PropertyName','PropertyValue'
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'plotbounds': logical indicating whether to plot bounds or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'maskEl': string indicating which element should be used to create
%   the preliminary standard mask [DEFAUL = 'C_13']
%   • 'ElLimits': nel x 2 matrix with rows corresponding to the
%   [min(el),max(el)] values for which the ylimits of the calibration
%   curves should be [DEFAULT = []]
%
% OUTPUTS ---------------------------------------------------------------
% Calstruct: Stucture with relevant output data
%   • Elnames:
%   • LineSlopeInt_Means: nel x 2 matrix with rows corresponding to the
%   [slope,intercept] values for the element corresponding to that row
%   • LineFitR2_Means = rsqmean;
%   • LineSlopeInt = p;
%   • LineFitR2 = rsq;
% dcellOut:
%   • d_div: L x 1 cell of structures where each structure is the
%        basic MIMS structure associated with the label
%   • line_types
%   • type_rows
%   • CalMasks
%   • dinMasks
%   • ppmVal
%   • CalLabels
%   • dmeans
%
%  Date           Author              E-mail                      Version
%  24 Mar  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
%   1 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%     Separated portions of the script into multiple functions since each
%     individually is commonly used

compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 1;
end

if strmatch('plotbounds',PropertyNames)
  bnds = PropertyVal{strmatch('plotbounds',PropertyNames)};
else
  bnds = 0;
end

if strmatch('saveimgs',PropertyNames)
  saveimgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  saveimgs = 0;
end

if strmatch('plotel',PropertyNames)
  plotel = interp_data_elemrng(d_div{1},PropertyVal{strmatch('plotel',PropertyNames)});
else
  plotel = interp_data_elemrng(d_div{1},[]);;
end

if strmatch('ellim',PropertyNames)
  ellims = PropertyVal{strmatch('ellim',PropertyNames)};
else
  ellims = [];
  %   ellims = [4e5,2e6;...'C_13'
  %     0e5,10e5;...'Ce140'
  %     0,1e6;...'Cu63'
  %     0,2e5;...'Cu65'
  %     0,6e5;...'Eu151'
  %     0,1e4;...'Fe57'
  %     0,1.5e5;...'Gd157'
  %     2.5e5,7.5e5;...'K_39'
  %     0,7e5;...'La138'
  %     0,3e5;...'Mn55'
  %     0,6e4;...'Mo95'
  %     1e5,7e5;...'P_31'
  %     0.5e4,2.5e4;...'S_33'
  %     2e4,9e4;...'S_34'
  %     3.5e7,5.5e7;...'S_36'
  %     0,5e5;...'Sr88'
  %     0,2e5;...'Ti48'
  %     0,6e5;...'Y_89'
  %     2e4,7e4;...'Zn64'
  %     0,2.25e4;...'Zn66'
  %   ];
end

%% Reads raw data
if isstr(d_div)
  [~,~,d_div,line_types,~,~,laONmask] = labelMIMSstruct(d_div);
end

delonly = rmfield(d_div{1},skip_fields);
elnames = sort(fieldnames( delonly ));
nel = length(elnames);

%% Creates a calibration mask using indicated channel
calibcells = find(cellfun(@(x) ~isempty(x),strfind(lower(line_types),'calib')));
nstnds = length(calibcells);
ppmnum = nan(nstnds,1);
caldata = cell(1,nstnds);
for c = 1:nstnds
  calmask = calibcells(c);
  dcal = d_div{calibcells(c)};
  dcalel = struct2cell(rmfield(dcal,skip_fields));
  calnumstr = strsplit(dcal.dataset,'_');
  ppmnum(c) = str2num(calnumstr{end});
  
  npnts = sum(laONmasks{calmask}(:));
  if npnts
    caldata{c} = reshape(cell2mat(cellfun(@(x) x(laONmasks{calmask}),dcalel,'uniformoutput',0)),...
      npnts,nel);
    
    %if plotimgs --- box plots
    %  figure('position',[284, 530, 1070, 427]);
    %  boxplot(dsmask{c}./repmat(max(dsmask{c}),npnts,1),elnames); hold on;
    %  hold on; plot(mean(dsmask{c})./max(dsmask{c}),'gx');
    %  title(sprintf('%s: Normalized box-whisker & means of data in Standard %d',strrep(dcal.dataset,'_','\_'),ppmnum(c)));
    %end
  else
    caldata{c} = nan(3,nel);
    disp(sprintf('%s does not have any included data in the mask',line_types{calmask}))
  end
end
[ppmnum,reorderrows] = sort(ppmnum);
caldata = caldata(reorderrows);
dmeans = cell2mat(...
  cellfun(@(x) nanmean(x)',caldata,'uniformoutput',0))';

dmedians = reshape(cell2mat(...
  cellfun(@(x) quantile(x,0.5)',caldata,'uniformoutput',0)),...
  nstnds,nel);

dstd = cell2mat(cellfun(@(x) nanstd(x)',caldata,'uniformoutput',0))';

dcellOut.dinMasks = caldata;
dcellOut.ppmVal = ppmnum;
dcellOut.CalLabels = line_types(calibcells(reorderrows));
dcellOut.dmeans = dmeans;
dcellOut.dstd = dstd;

%% Linear curve fitting
dsmaskmat = cell2mat(caldata(:));
nonNAN = ~isnan(dsmaskmat);
nonNANmeans = ~isnan(dmeans);
x = cell(nstnds,1);
for c = 1:nstnds
  x{c} = repmat(ppmnum(c),size(caldata{c},1),1);
end
x = cell2mat(x);

p = nan(nel,2);
rsq = nan(nel,1);
%S = cell(1,nel);
%mu = cell(1,nel);
for i = 1:nel
  [fobj,gof,outstruct] = fit(x(nonNAN(:,i)),dsmaskmat(nonNAN(:,i),i),'poly1');%'Lower', [-Inf dmeans(1,i)], 'Upper', [Inf dmeans(1,i)]);
  p(i,:) = [fobj.p1,fobj.p2];
  rsq(i) = gof.rsquare;
  %[p(i,:),S{i},mu{i}] = polyfit(x,dsmaskmat(:,i),1);
end

pmeans = nan(nel,2);
rsqmean = nan(nel,1);
for i = 1:nel
  [fobj,gof,outstruct] = fit(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),'poly1');%'Lower', [-Inf dmeans(1,i)], 'Upper', [Inf dmeans(1,i)]);
  pmeans(i,:) = [fobj.p1,fobj.p2];
  rsqmean(i) = gof.rsquare;
end
datastr = strtok(d_div{calibcells(1)}.dataset,'Cal');
datastr = datastr(1:end-1);

%% Depicting linear fit with raw data points
if plotimgs
  for i = plotel(:)'%1:nel
    ccel = figure('name',elnames{i});
    plot(x,dsmaskmat(:,i),'.c');
    hold on;
    if bnds
      boundedline(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),dstd(nonNANmeans(:,i),i),'b', 'alpha');
      boundedline(ppmnum(nonNANmeans(:,i)),dmeans(nonNANmeans(:,i),i),dstd(nonNANmeans(:,i),i),'xr', 'alpha');
    else
      plot(ppmnum,p(i,1).*ppmnum + p(i,2),'b');
      plot(ppmnum,dmeans(:,i),'xr');
      plot(ppmnum,pmeans(i,1).*ppmnum + pmeans(i,2),'r');
    end
    
    title(sprintf('%s: %s Calibration curves',...
      strrep(datastr,'_','\_'),strrep(elnames{i},'_','\_')));
    if isempty(ellims)
      mxyval = quantile(dsmaskmat(:,i),0.98);
      ylim([min(dsmaskmat(:,i)),mxyval]);
    else
      mxyval = ellims(i,2);
      ylim(ellims(i,:));
    end
    text(0,mxyval*.95,...
      sprintf('$$I = %1.2f C + %1.2f, R^2 = %1.4f$$',p(i,1),p(i,2),rsq(i)),...
      'fontsize',14,'color','b','VerticalAlignment','top','interpreter','latex')
    text(0,mxyval*.9,...
      sprintf('$$I = %1.2f C + %1.2f, R^2 = %1.4f$$',pmeans(i,1),pmeans(i,2),rsqmean(i)),...
      'fontsize',14,'color','r','VerticalAlignment','top','interpreter','latex')
    xlabel('C = Standard Concentration [ppm]')
    ylabel('I = Intensity [cps]')
  end
  if saveimgs
    close all
    save_open_figures(strcat(compdir.lafldr,datastr),[],[],'CalibrationCurve')
    %export_fig(strcat('Images\CalibrationCurve_',elnames{i}),...
    %  '-png');
    %close(ccel)
  end
end

dockf on all
Calstruct.Elnames = elnames;
Calstruct.LineSlopeInt_Means = pmeans;
Calstruct.LineFitR2_Means = rsqmean;
Calstruct.LineSlopeInt = p;
Calstruct.LineFitR2 = rsq;

end
