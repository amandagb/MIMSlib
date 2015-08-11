function [la0info,varargout] = find_laoff_ind(data,elem_ana_range,varargin);

%% Script: [laoff_ind,varargout] = find_laoff_ind(data,elem_ana_range,varargin);
% Description:
% Example:  ila0 = find_laoff_ind(dA,[],'ana',0);
%           [ila0,dAla0] = find_laoff_ind(dA,[],'ana',0);
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'ana': plots line by line data of laser off region and gives mean, std,
%   and rsd
% OUTPUTS ---------------------------------------------------------------
% la0info: structure which contains information about laser off periods
%   la0ind:
%   meanla0ind:
%   la0peak:
%   meanla0peak:
%
%  Date           Author            E-mail                      Version
%  25 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  30 Aug 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  13 Sept 2011  Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%  18 June 2012  Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Tries to match the resulting image to either a step function
%     (forward or backward) or a rectangle. Signal goes from low-high,
%     high-low, or low-high-low
%  19 June 2012  Amanda Gaudreau   amanda.gaudreau@gmail.com     4.1
%     Uses information from mean dataset to guide numbers for overall
%     datasets
%  24 July 2012  Amanda Gaudreau   amanda.gaudreau@gmail.com     4.2
%     Replaced subfunctions with find_change_thresh function which was
%     created as an external function for finding changes in vector data

regions = [];
fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
allflds = fieldnames(data);
rmflds = allflds(ismember(allflds,fldnm)==0);
nel = length(fldnm); 

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('ana',PropertyNames)
  ana = PropertyVal{strmatch('ana',PropertyNames)};
else
  ana = 0;
end

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end

% Removing nans is required for 120710Brain dataset as of 9/16/2011 -- not
% incorporated in newest version (v5)
if isstruct(data)
  d = structfun(@(x) normalize_data(x,[min(x(:)),max(x(:))]),rmfield(data,rmflds),'uniformoutput',0);
  [M,N] = size(getfield(data,fldnm{1}));
  d = reshape(struct2array(d),[M,N,nel]);
else 
  d = normalize_data(data,[min(data(:)),max(data(:))]);
end
regions = setfield(regions,'data',d);
data2D = sum(d,3);
[M,N] = size(data2D);

% Find overall numbers
meandata1D = mean(data2D,1); % comment out to find line by line information
smoothmean = medfilt1(meandata1D,9);
sgndiff_mean = sign(diff(smoothmean));
[meanla0ind,meanchanges,meantopthresh] = find_change_thresh(meandata1D,'starting_point','t'); 

meanla0peak = zeros(1,2);
if meanla0ind(1)
  front_ind = 1:round(meanla0ind(1)*1.5);
  m1sngs = find(sgndiff_mean(meanla0ind(1):end) ~= 1) + meanla0ind(1);
  meanla0peak(1) = m1sngs(1) - 1;
else
  front_ind = [];
end

if meanla0ind(2)
  backdiff = round((N-meanla0ind(end))*1.5);
  backstart = N - backdiff;
  back_ind = backstart:N;  
  flr_sgndiff_datavec =  fliplr(sgndiff_mean);
  equiv_ind = length(meandata1D)-meanla0ind(2)+ 1;
  m1sngs = find(flr_sgndiff_datavec(equiv_ind:end) ~= -1) + equiv_ind;
  meanla0peak(2) = length(meandata1D)-m1sngs(1)+1+1;
else
  back_ind = [];
end

% Do same analysis for each line using information from the mean data
% vector
top_ind = zeros(M,2);
bottom_ind = zeros(M,2);
endstartrise = zeros(M,1);
beginendfall = zeros(M,1);
mintopchanges = ones(M,1).*100;
minbottomchanges = ones(M,1).*100;
if p; h = figure; end
for k = 1:size(data2D,1)
  datavec = data2D(k,:);
  smoothdatavec = medfilt1(datavec,9);
  sgndiff_datavec = sign(diff(smoothdatavec));
  shortdata = datavec([front_ind,back_ind]);
  Ns = length(shortdata);
  [top_ind(k,:),mintopchanges(k),intermtopthresh] = find_change_thresh(shortdata,'starting_point','t'); 
  top_ind(k,2) = N - (Ns - top_ind(k,2));
  
  [bottom_ind(k,:),minbottomchanges(k),intermbottomthresh] = find_change_thresh(shortdata,'starting_point','b');
  if bottom_ind(k,2); bottom_ind(k,2) = N - (Ns - bottom_ind(k,2)); end
  
  if top_ind(k,1) ~= 0
    m1sngs = find(sgndiff_datavec(top_ind(k,1):end) ~= 1) + top_ind(k,1);
    endstartrise(k) = m1sngs(1) - 1;
  end
  
  if top_ind(k,2) ~= 0
    flr_sgndiff_d_comb_vec =  fliplr(sgndiff_datavec);
    equiv_ind = length(datavec)-top_ind(k,2)+ 1;
    m1sngs = find(flr_sgndiff_d_comb_vec(equiv_ind:end) ~= -1) + equiv_ind;
    beginendfall(k) = length(datavec)-m1sngs(1)+1+1;
  end
  
  if p
    plot(datavec,'.'); hold on;
    plot(top_ind(k,:)',datavec(top_ind(k,:))','go');
    plot(endstartrise(k),datavec(endstartrise(k)),'ro');
    plot(beginendfall(k),datavec(beginendfall(k)),'ro');
    title(sprintf('Line %d',k));
    hold off;
  end
end
la0peak = [endstartrise,beginendfall];
la0ind = top_ind;

% Construct output variable
la0info.top_la0ind = top_ind;
la0info.bottom_la0ind = bottom_ind;
la0info.meanla0ind = meanla0ind;
la0info.la0peak = la0peak;
la0info.meanla0peak = meanla0peak;

if ana
  la0dat = [];
  info = [];
  fldnm = fieldnames(data);
  for f = 1:length(fldnm)
    if isstruct(data)
      d = getfield(data,fldnm{f});
    else
      d = data;
    end
    col = jet(size(d,1));
    if length(fldnm{f}) == 6 %isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
      info = setfield(info,'data',d(:,i_la0));
      info = setfield(info,'mean',mean(info.data(:)));
      info = setfield(info,'std',std(info.data(:)));
      info = setfield(info,'rsd',info.std/info.mean);
      
      la0dat = setfield(la0dat,fldnm{f},info);
      if p
        figure('color','w');
        hold on;
        lnstr = cell(0,size(d,1));
        for i = 1:size(d,1)
          plot(info.data(i,:),'color',col(i,:),'linewidth',2,'displayname',sprintf('line %1g',i));
          lnstr{i} = sprintf('line %1g',i);
        end
        %legend(lnstr);
        title(sprintf('%s LA off data summary, samples %d - %d',fldnm{f},la0ind(1),la0ind(end)),'fontsize',fntsz);
        xlabel('Sample #','fontsize',fntsz);
        ylabel('Intensity (cps)','fontsize',fntsz);
        x_rng = get(gca,'xlim');
        y_rng = get(gca,'ylim');
        text(x_rng(end) - 3,((5*(y_rng(2)-diff(y_rng)))/6)+y_rng(1),...
          sprintf('mean = %d \n \sigma = %d \n rsd = %d',info.mean,info.std,info.rsd),'fontsize',fntsz);
      end
    end
  end
  varargout{1} = la0dat;
  dockf on all;
else varargout{1} = [];
end

%% Subfunctions
  function [thresh_ind,min_changes,interm_thresh,la0ind] = findudrise(dvec);
    Ddvec = diff(dvec);
    sortdvec = sort(dvec);
    min_changes = 100;
    thresh_ind = [0,0];
    interm_thresh = 0;
    
    top_thresh = mean(dvec);
    [v,top_threshi] = min(abs(sortdvec-top_thresh));
    top_thresh = sortdvec(top_threshi);
    highlowfail = 1;
    top10vec = zeros(size(dvec));
    while highlowfail && top_threshi > 1
      top10vec = dvec > top_thresh;
      chngs = find(diff(top10vec) ~= 0);
      if length(chngs) == 1
        highlowfail = 0;
        la0ind = find(top10vec == 0);
        if top10vec(chngs(1)) == 0
          thresh_ind = [chngs,0];
        else
          thresh_ind = [0,chngs];
        end
      elseif length(chngs) == 2 && top10vec(chngs(1)) == 0
        highlowfail = 0;
        la0ind = find(top10vec == 0);
        thresh_ind = chngs;
      else
        top_threshi = top_threshi - 1;
        top_thresh = sortdvec(top_threshi);
      end
      
      if length(chngs) < min_changes && rem(length(chngs),2) == 0
        min_changes = length(chngs);
        interm_thresh = sortdvec(top_threshi);
      end
    end
    if rem(length(chngs),2) == 0;
      min_changes = min(length(chngs),min_changes);
      interm_thresh = top_thresh;
    end
  end

  function [thresh_ind,min_changes,interm_thresh,la0ind] = finddurise(dvec);
    Ddvec = diff(dvec);
    sortdvec = sort(dvec);
    min_changes = 100;
    thresh_ind = [0,0];
    interm_thresh = 0;
    
    bottom_thresh = sortdvec(4);
    bottom_threshi = 4;
    bottom10vec = zeros(size(dvec));
    highlowfail = 1;
    while highlowfail && bottom_threshi < length(dvec)
      bottom10vec = dvec > bottom_thresh;
      chngs = find(diff(bottom10vec) ~= 0);
      if length(chngs) == 1
        highlowfail = 0;
        la0ind = find(bottom10vec == 0);
        if bottom10vec(chngs(1)) == 0
          thresh_ind = [chngs,0];
        else
          thresh_ind = [0,chngs];
        end
      elseif length(chngs) == 2 && bottom10vec(chngs(1)) == 0
        highlowfail = 0;
        la0ind = find(bottom10vec == 0);
        thresh_ind = chngs;
      else
        bottom_threshi = bottom_threshi + 1;
        bottom_thresh = sortdvec(bottom_threshi);
      end
      
      if length(chngs) < min_changes && rem(length(chngs),2) == 0 ...
          || length(chngs) == 1
        min_changes = length(chngs);
        interm_thresh = sortdvec(bottom_threshi);
      end
      
    end
    if rem(length(chngs),2) == 0;
      min_changes = min(length(chngs),min_changes);
      interm_thresh = bottom_thresh;
    end
  end

end