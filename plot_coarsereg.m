function outputs = plot_coarsereg(metric_str,phi,param_cell,param_names,varargin)
%% Script:  
% Description:  
% Example:  plot_coarsereg('MI',phi_MI,param_cell,param_names,'normalize',1);
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'color_res':  numeric specifying color resolution of the heatmap
%   (maximum 256) [DEFAULT = 256]
%   'thresh_data':  1 x 2 matrix with minimum threshold value as first
%   element and maximum threshold value as second
%   'lines': matrix containing line numbers user would like to plot
%   [DEFAULT: ALL LINES]
%   'samples': matrix containing samples (S) user would like to plot [DEFAULT:
%   ALL SAMPLES]
%   'plot_type': string - 'log','linear','both' [DEFUALT: BOTH]
%   'elim_blanks': 1 or 0 - eliminates any lines which have NAN data
% OUTPUTS ---------------------------------------------------------------
% 
%
%  Date           Author            E-mail                      Version
%  14 Aug  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('normalize',PropertyNames)
  norm_met = PropertyVal{strmatch('normalize',PropertyNames)};
else
  norm_met = 0;
end

xtran = reshape(param_cell{1}',1,size(param_cell{1},1)*size(param_cell{1},2)); 
nxtran = size(param_cell{1},2);
[sxtran,sxtranid] = sort(xtran);
ytran = reshape(param_cell{2}',1,size(param_cell{2},1)*size(param_cell{2},2)); 
nytran = size(param_cell{2},2);
[sytran,sytranid] = sort(ytran);
deg_ang = reshape(param_cell{3}',1,size(param_cell{3},1)*size(param_cell{3},2)); 
ndeg_ang = size(param_cell{3},2);
[sdeg,sdegid] = sort(deg_ang);
xsc = reshape(param_cell{4}',1,size(param_cell{4},1)*size(param_cell{4},2));   
nxsc = size(param_cell{4},2);
[sxsc,sxscid] = sort(xsc);
ysc = reshape(param_cell{5}',1,size(param_cell{5},1)*size(param_cell{5},2));  
nysc = size(param_cell{5},2);
[sysc,syscid] = sort(ysc);
sk = reshape(param_cell{6}',1,size(param_cell{6},1)*size(param_cell{6},2)); 
nsk = size(param_cell{6},2);
[ssk,sskid] = sort(sk);
%sparamid = [sxtranid;sytranid;sdegid;sxscid;syscid;sskid];
row_width = nxsc*nysc*nsk;%size(param_cell{4},2)*size(param_cell{5},2)*size(param_cell{6},2);
col_width = nxtran*nytran*ndeg_ang;%size(param_cell{1},2)*size(param_cell{2},2)*size(param_cell{3},2);
nruns = size(param_cell{1},1);
param_length = [nxtran,nytran,ndeg_ang,nxsc,nysc,nsk];

ind2paramvec = zeros(row_width*col_width,6);
niter = 0;
for a = 1:nxtran %transx (25's - H)
  for b = 1:nytran %transy (5's - H)
    for c = 1:ndeg_ang %rot (1's - H)
      for d = 1:nxsc %scalex (25's - V)
        for e = 1:nysc %scaley (5's - V)
          for f = 1:nsk %skew (1's - V)
            niter = niter + 1;
            ind2paramvec(niter,:) = [a,b,c,d,e,f];
          end
        end
      end
    end
  end
end

transparamsind = zeros(nruns,6);
transparams = zeros(nruns,6);
dep_varparam = zeros(6,size(param_cell{1},1)*size(param_cell{1},2));
lncolrs = colormap(lines(6));
for c = 1:nruns
  figure; 
  imagesc(reshape(phi(:,c),row_width,col_width)); axis image;
  title(sprintf('%s: Grid %d',metric_str,c)); 
  hold on;
  for j = 1:size(xsc,2)-1 
    plot([0,col_width],j.*size(ysc,2)*size(sk,2).*ones(1,2)+0.5,'--w');
  end
  
  for j = 1:size(xtran,2)-1
    plot(j.*size(ytran,2)*size(deg_ang,2).*ones(1,2)+0.5,[0,row_width],'--w');
  end
  
  [maxval,maxind] = max(phi(:,c));
  xind = ceil(maxind/row_width);
  yind = maxind - (xind-1)*row_width;
  hold on; plot(xind,yind,'.w');
  
  scalexind = ceil(yind/(nysc*nsk));
  scaleyind = ceil((yind - nysc*nsk*(scalexind-1))/nsk);
  skewind = rem(yind,nsk) + nsk*0^rem(yind,nsk);
  transxind = ceil(xind/(nytran*ndeg_ang));
  transyind = ceil((xind - (nytran*ndeg_ang)*(transxind-1))/ndeg_ang);
  rotind = rem(xind,ndeg_ang) + ndeg_ang*0^rem(xind,ndeg_ang);
  transparamsind(c,:) = [transxind, transyind, rotind, scalexind, scaleyind, skewind];
  transparams(c,:) = [xtran(nxtran*(c-1) + transxind), ...
    ytran(nytran*(c-1) + transyind), ...
    deg_ang(ndeg_ang*(c-1) + rotind),...
    xsc(nxsc*(c-1) + scalexind),... 
    ysc(nysc*(c-1) + scaleyind),...
    sk(nsk*(c-1) + skewind)];
  
  for k = 1:6
    ind_varparamind = transparamsind(c,:);
    for i = 1:param_length(k)
      findparams = ind_varparamind;
      findparams(k) = i;
      dep_varparam(k,i + param_length(k)*(c-1)) = ...
        phi(find(sum(abs(ind2paramvec - repmat(findparams,size(ind2paramvec,1),1)),2) == 0),c);
    end
    if norm_met
      dep_varparam(k,[1:param_length(k)] + param_length(k)*(c-1)) = ...
        dep_varparam(k,[1:param_length(k)] + param_length(k)*(c-1))./...
        max(dep_varparam(k,[1:param_length(k)] + param_length(k)*(c-1)));
    end
  end
end

for i = 1:6
  figure;
  %plot(param_cell{i}',reshape(dep_varparam(i,1:param_length(i))',param_length(i),nruns),'marker','.');
  plot(param_cell{i}',reshape(dep_varparam(i,:)',param_length(i),nruns),'marker','.');
  hold on;
  for c = 1:nruns
    hobj = plot(param_cell{i}(c,transparamsind(c,i)),...
      dep_varparam(i,transparamsind(c,i) + param_length(i)*(c-1)),...
      'o','markerfacecolor',lncolrs(c,:),...
      'markeredgecolor',lncolrs(c,:));
    hAnnotation = get(hobj,'Annotation');
    hLegendEntry = get(hAnnotation,'LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
  end
  title(sprintf('%s: %s',metric_str, param_names{i}));
  axis tight
  ylabel(metric_str); xlabel(param_names{i});
end

figure;
for i = 1:6
  plot(dep_varparam(i,1:param_length(i)),'color',lncolrs(i,:),'marker','.','displayname', param_names{i});
  hold on;
  for c = 1:nruns
    hobj = plot(transparamsind(c,i) + param_length(i)*(c-1),...
      dep_varparam(i,transparamsind(c,i) + param_length(i)*(c-1)),...
      'o','markerfacecolor',lncolrs(i,:),'markeredgecolor',lncolrs(i,:));
    hAnnotation = get(hobj,'Annotation');
    hLegendEntry = get(hAnnotation,'LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
  end
end
legend(param_names,'location','north')
title(sprintf('%s: all parameters & runs',metric_str));
dockf on all


