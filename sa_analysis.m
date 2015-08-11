function datasummary = sa_analysis(varargin)
%% Script:  d0314 = sa_analysis('datestr','03_15','nvar',3,'plotconv',1);
%           dsum = sa_analysis('datestr','05_09','nvar',6,'plottype','lin','plotdetail','m','Fimg',Fimg,'Mimg',Mimg,'Medgeth',5,'plotMI',1);
% Description: Reads data from MM_DD_nvar6Paccepth1_T#_{6varBinary/Image label}_#ROWS
% Given an input
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'datestr'
%   'nvar'
%   'plottype': string indicating whether to plot a 3-D plot ('3d' for 2 -
%     3 varaibles) or 2-D slice plots ('2d' for 3 variables). For > 3
%     parameters, use 'lin' [DEFAULT = '' (no plot)]
%   'plotdetail': string indicating whether to plot markers 'm'
%     corresponding to reannealing [DEFAULT don't use markers]
%   'plotcnvg': binary indicating whether to plot convergence plots for
%     each of the reannealing intervals in the dataset [DEFAULT = 0]
%   'plotMI': binary indiciating whether to plot MI distribution plots or
%   not [DEFAULT = 0]
%   'saveplot': binary indicating whether to save plots or not [DEFAULT = 0]
%   'Mimg': original moving image before resizing and before registration
%     parameters are applied
%   'Fimg': original fixed image before resizing
%     NOTE ON IMAGE: For the metallomic image, include the entire
%     structure. The appropriate field will be extracted
%   'Medgeth': Moving image threshold used for overlaying function
%
% OUTPUTS ---------------------------------------------------------------
% datasummary:  D x 3 cell where each row contains data from a separate
%     experiment. The columns have the following information:
%     COL 1: string containing the name of the dataset csv file analyzed
%     COL 2: Nre x  matrix with the following information where Nre is the
%     number of reannealing runs
%       col 1:nvar -- best state in reannealing interval
%       col nvar + 1 -- "best" (minimum) MI range in the interval
%       col nvar + 2 -- iteration of the minimum MI range
%       col nvar + 3 -- percentage of points in interval where the MI value was < some MI thresh
%       col nvar + 4 -- endpntsltMIthresh
%       col nvar + 5 -- iter2conv
%       col nvar + 6 -- total number of iterations in interval
%       col nvar + 7 -- temperature to convergence
%     COL 3: Structure with the following fields:
%       mu0, MI0, muexp, MIexp, T0, tau, naccept, mu_num, bounds,
%       scaled_bounds
%
%  Date           Author            E-mail                      Version
%  14 Mar  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  19 Mar  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Updating function for reading in multimodal data -- only difference
%     should be file names
%  15 May  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Cleaned up commented code and used newly written imoverlay.m function
%     for plotting overlay images. This function is multipurpose & much
%     more robust
%  21 Oct  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Changed the naming method of the data saved so the string has been
%     updated to reflect the new naming scheme
%  12 May  2015   Amanda Gaudreau   amanda.gaudreau@gmail.com     5

set(0,'defaultfigurecolor','w','DefaultAxesFontSize',14,...'DefaultAxesFontWeight','bold','DefaultAxesLineStyleOrder','-|--|:|-.',...
  'DefaultLineLineWidth',2,... 'Defaultaxesposition','remove',...
  'defaultAxesFontName','Arial');

cdir = loaddirfun;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('datestr',PropertyNames)
  mmddstr = PropertyVal{strmatch('datestr',PropertyNames)};
else
  mmddstr = '03_05';
end

if strmatch('nvar',PropertyNames)
  nvar = PropertyVal{strmatch('nvar',PropertyNames)};
else
  nvar = [];
end
nvarstr = strcat('nvar',num2str(nvar));

if strmatch('T0',PropertyNames)
  T0 = PropertyVal{strmatch('T0',PropertyNames)};
else
  T0 = [];
end
T0str = strcat('_T',num2str(T0));

usemrk = 0;
if strmatch('plotdetail',PropertyNames)
  plotdetail = PropertyVal{strmatch('plotdetail',PropertyNames)};
  if strmatch('m',plotdetail)
    usemrk = 1;
  end
else
  usemrk = 0;
end

if strmatch('plottype',PropertyNames)
  plottype = PropertyVal{strmatch('plottype',PropertyNames)};
else
  plottype = '';
end

if strmatch('plotcnvg',PropertyNames) % plot reannealing convergence
  pana = PropertyVal{strmatch('plotcnvg',PropertyNames)};
else
  pana = 0;
end

if strmatch('plotMI',PropertyNames) % plot reannealing convergence
  pMI = PropertyVal{strmatch('plotMI',PropertyNames)};
else
  pMI = 0;
end

if strmatch('plotFnum',PropertyNames) % plot reannealing convergence
  pF = PropertyVal{strmatch('plotFnum',PropertyNames)};
else
  pF = [];
end

if strmatch('saveplot',PropertyNames)
  svp = PropertyVal{strmatch('saveplot',PropertyNames)};
else
  svp = 0;
end

if strmatch('Mimg',PropertyNames)
  Mimg = PropertyVal{strmatch('Mimg',PropertyNames)};
else
  Mimg = [];
end

if strmatch('Fimg',PropertyNames)
  Fimg = PropertyVal{strmatch('Fimg',PropertyNames)};
else
  Fimg = [];
end

mu0 = [0,0,0,1,1,0];
currentfldr = cd;
Dpos = strfind(currentfldr,'Documents');
comp = currentfldr(strfind(currentfldr,'Users\')+6:Dpos-2);%'TENN1_000'; % MED comp %'ADGB'; % ThinkPad; %
dboxpath = strcat(currentfldr(1:Dpos-1),'Documents\Dropbox');
path = strcat(dboxpath,'\MADLab Research\Data\simanneal Experiments\');%'D:\My Documents\Dropbox\MADLab Research\Data\simanneal Experiments\';%
paramstr = {'tx','ty','\theta','sx','sy','sk'};
markers = {'o','s','d','^','h'};
ncolbins = 16;
nslices = 15;
varbinstr = '';
pctMIthresh = 5;
finalpctofiter = 10;
pntthreshltMIthresh = 90;
mmddslashi = length(mmddstr)+1;

filestrs = {'Paccepth1','genfctn','besth1','png'};
F = ls(path);
Fcell = cellstr(F);
% Isolates slashes and csv files
Fslashi = strfind(Fcell,'_');
Fcsvi = strfind(Fcell,'.csv');
Fanai = strfind(Fcell,filestrs{1});
relFi = cellfun(@(x) isempty(x) == 0,Fslashi);
relFii = cellfun(@(x) isempty(x) == 0,Fcsvi);
relFiii = cellfun(@(x) isempty(x) == 0,Fanai);
relF = find(relFi + relFii + relFiii== 3);
F = F(relF,:);
Fcell = Fcell(relF);
Fslashi = Fslashi(relF);
Fcsvi = Fcsvi(relF);

% Isolate strings which indicate files to use
if ~isempty(mmddstr)
  relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,mmddstr));
else
  relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,'csv'));
end
relFii = cellfun(@(x) isempty(x) == 0,strfind(Fcell,nvarstr));
relFiii = cellfun(@(x) isempty(x) == 0,strfind(Fcell,T0str));
if ~isempty(varbinstr)
  relFiv = cellfun(@(x) isempty(x) == 0,strfind(Fcell,varbinstr));
else
  relFiv = cellfun(@(x) isempty(x) == 0,strfind(Fcell,'csv'));
end
relF = find(relFi + relFii + relFiii + relFiv == 4);
F = F(relF,:);
Fcell = Fcell(relF);
Fslashi = Fslashi(relF);
Fcsvi = Fcsvi(relF);
MM_DDarry = F(:,1:mmddslashi-1);
MM_DD = cellstr(unique(MM_DDarry,'rows'));

tempi = cellfun(@(x) x + 2,strfind(Fcell,'_T'));
nvari = cellfun(@(x) x + 4,strfind(Fcell,'nvar'));
varbini = cellfun(@(x) x(3)+1,Fslashi);%cellfun(@(x) x - 6,strfind(Fcell,'.csv'));

lenF = length(Fcell);
varbinmat = zeros(lenF,6);
imrows = zeros(lenF,1);
MIMStag = zeros(lenF,1);
elused = cell(lenF,1);
S = cell(1,lenF);
S(:) = {'                              '};
varbinchar = char(S);
% varbinchar = cell(lenF,1);
for f = 1:lenF
  if varbini(f) - 1 == Fslashi{f}(end)
    varbinchar(f,:) = F(f,varbini(f):(Fcsvi{f}-1));
    imrows(f) = 32;
  else
    str = F(f,varbini(f):(Fslashi{f}(end)-1));
    varbinchar(f,1:length(str)) = str;
    imrows(f) = str2num(strrep(F(f,(Fslashi{f}(end)+1):(Fcsvi{f}-1)),'i',''));
  end
  if length(varbinchar) == 6
    activeparams = strfind(varbinchar(f,:),'1');
    varbinmat(f,activeparams) = 1;
  else
    varbinmat(f,:) = 1;
    MIMStag(f) = 0;%str2num(F(f,Fslashi{f}(end-1)-1));%(Fslashi{f}(4)+1):(Fslashi{f}(5)-1)));
    elused{f} = F(f,(Fslashi{f}(end-1)+1):(Fslashi{f}(end)-1));  varbinchar(f,strfind(varbinchar(f,:),'_')+1:end);
  end
end
unqvarbin = unique(varbinchar,'rows');
reorderi = [];
Ftemp = [];
for u = 1:size(unqvarbin,1)
  matchvarbin = sum(abs(varbinchar - repmat(unqvarbin(u,:),lenF,1)),2);
  rei = find(matchvarbin == 0);
  for t = 1:length(rei)
    reli = rei(t);
    eqslashi = find(Fslashi{reli} == tempi(reli)-2);
    tempstr = F(reli,tempi(reli):(Fslashi{reli}(eqslashi+1)-1));
    tempall(t) = str2num(regexprep(regexprep(tempstr,'p','.'),'i',''));
  end
  [tempunq,reordert] = sort(tempall);
  reorderi = [reorderi;rei(reordert)];
  Ftemp = [Ftemp(:);tempunq(:)];
  clear tempall
end
Fcell = Fcell(reorderi);
F = F(reorderi,:);
c = sastatescalefactors;
imrows = imrows(reorderi);
MIMStag = MIMStag(reorderi);;
elused = elused(reorderi);
datasummary = cell(lenF,3);
datasummary(:,1) = Fcell;
for f = 1:lenF
  fname = Fcell{f};
  fstr = fname;
  Fmmdd = MM_DDarry(f,:);
  % If date was before Mar 6, use the following indices
  if str2num(Fmmdd(1:2)) < 3
    MIpasti = 1; MIpropi = 2; statepropi = 9:(8+nvar); statepasti = 9:(8+nvar);
    delEi = 3; h1i = 4; Ui = 6; Ti = 7; acci = 8;
  elseif str2num(Fmmdd(4:end)) < 6 && str2num(Fmmdd(1:2)) == 3
    MIpasti = 1; MIpropi = 2; statepropi = 9:(8+nvar); statepasti = 9:(8+nvar);
    delEi = 3; h1i = 4; Ui = 6; Ti = 7; acci = 8;
  else
    % If date was on or after Mar 6, use the following indices
    MIpasti = 1; statepasti = [2:nvar+1]; MIpropi = nvar + 2; statepropi = [(nvar + 3):(2*nvar + 2)];
    delEi = 2*nvar + 3; h1i = 2*nvar + 4; Ui = 2*nvar + 5; Ti = 2*nvar + 6; acci = 2*nvar + 7;
  end
  
  slashi = strfind(fname,'_');
  if isempty(elused{f})
    binarycode = fname(slashi(4)+1:slashi(4)+6);
  else
    binarycode = '111111';%ones(1,6);
  end
  fnameheader = fname(1:slashi(2)+5); % MM_DD_nvar#
  fnamefooter = fname(slashi(4)+1:end-4); % variable binary code
  for j = 1:length(slashi)
    fstr = strcat(fstr(1:(slashi(j)-1 + (j-1))),'\',fstr((slashi(j) + (j-1)):end));
  end
  activeparams = strfind(binarycode,'1');
  d = importdata(strcat(path,fname));
  if any(isnan(d(end,:)))
    Fm = d(end,1); Fn = d(end,2);
    time_hr = d(end,3)/3600;
    d = d(1:end-1,:);
  else
    Fm = 32; Fn = 32; time_hr = nan;
  end
  ub = zeros(1,length(activeparams)); lb = zeros(1,length(activeparams));
  ub(activeparams) = max(d(:,statepropi));%max[Fn*0.2,Fm*0.2,10*pi/180,1.3,1.3,0.2].*c;%ub = [10,10,10*pi/180,1.3,1.3,0.2].*c;%
  lb(activeparams) = min(d(:,statepropi));
  c = c(activeparams);
  pstr = regexprep(sprintf('%s_T%d_%s',fnameheader,Ftemp(f),fnamefooter),'_','\\_');
  
  % Reading out variables from dataset
  initstate = d(1,statepasti);
  MIpast = d(:,MIpasti);
  MIprop = d(:,MIpropi);
  accvec = logical(d(:,acci));
  MIused = MIpast;
  MIused(accvec) = MIprop(accvec);
  stateused = d(:,statepasti);
  stateused(accvec,:) = d(accvec,statepropi);
  npnts = sum(accvec);
  accstate = d(accvec,statepropi);
  accstate = accstate./repmat(c,size(accstate,1),1);
  accstatei = find(accvec == 1);
  MIstate = abs(MIprop(accvec));
  Tvec = d(:,Ti);
  dTvec = diff(Tvec);
  T0 = Tvec(1);
  sgnTvec = sign(dTvec);
  reannealstart = [1;find(sgnTvec == 1) + 1;length(Tvec)];
  annpnts = diff(reannealstart);
  increannind = find(annpnts > 1000)+1;
  reannealstart = reannealstart([1,increannind']);
  nreanneals = length(reannealstart) - 1;
  if length(reannealstart) < 2
    tempdecay = mean(Tvec(2:21)./Tvec(1:20));
  else
    tempdecay = mean(Tvec(2:reannealstart(2)-1)./Tvec(1:reannealstart(2)-2));
    naccept = sum(accvec(1:reannealstart(2)));
  end
  % Evaulation Parameters
  maxMI = max(MIused);
  [bestMI,besti] = min(MIused);
  beststate = d(besti,statepropi);
  MIrng = bestMI - maxMI;
  pctMI = bestMI - pctMIthresh*MIrng/100;
  iterltpct = find(MIused < pctMI);
  nannealpnts = zeros(nreanneals,1);
  beststateinrng = zeros(nreanneals,nvar);
  MIrngmin = zeros(nreanneals,1);
  bestinrngk = zeros(nreanneals,1);
  totpntltMIthresh = zeros(nreanneals,1);
  endpntsltMIthresh = zeros(nreanneals,1);
  iter2conv = zeros(nreanneals,1);
  temp2conv = zeros(nreanneals,1);
  for r = 1:nreanneals
    nannealpnts(r) = reannealstart(r+1)-reannealstart(r);
    beginendbuffer = round(nannealpnts(r)*finalpctofiter/100);
    reannealind = reannealstart(r):(reannealstart(r+1)-1);
    relind = (reannealstart(r)+beginendbuffer):(reannealstart(r+1)-1);
    % Best MI value & corresponding state for each reannealing interval
    [MIrngmin(r),bestinrngi] = min(MIused(relind));
    bestinrngk(r) = bestinrngi;
    bestinrngi = relind(bestinrngi);
    beststateinrng(r,:) = stateused(bestinrngi,:);
    % Look for relative indices of points below the MI threshold
    relpcti = find(iterltpct >= reannealstart(r) &...
      iterltpct < reannealstart(r+1));
    relpct = iterltpct(relpcti);
    totpntltMIthresh(r) = length(relpct)/nannealpnts(r);
    pctind = reannealstart(r+1) - beginendbuffer - 1;
    endpntsltMIthresh(r) = sum(relpct > pctind)/beginendbuffer;
    indltthreshlog = ismember(reannealind,relpct);
    sumindlt = fliplr(cumsum(fliplr(indltthreshlog))./[1:nannealpnts(r)]);
    indltpntthresh = find(sumindlt >= pntthreshltMIthresh/100);
    if isempty(indltpntthresh)
      iter2conv(r) = nan;
      temp2conv(r) = nan;
    elseif indltpntthresh(1) == nannealpnts(r)
      iter2conv(r) = nan;
      temp2conv(r) = nan;
    else
      iter2conv(r) = indltpntthresh(1) - 1;
      temp2conv(r) = Tvec(reannealind(indltpntthresh(1)));%
    end
    if pana
      figure; plot(MIused(reannealind));
      hold on; plot([0,length(reannealind)],pctMI.*ones(1,2),':k');
      plot((nannealpnts(r)-beginendbuffer).*ones(1,2),[maxMI,bestMI],':k');
      plot((iter2conv(r)).*ones(1,2),[maxMI,bestMI],':r');
      title(sprintf('%s MI',pstr));
      ylabel(sprintf('MI (range: %1.2f - %1.2f)',maxMI,bestMI));
    end
  end % END reannealing loop
  muinit = mu0; muinit(activeparams) = initstate./c;
  mubest = mu0; mubest(activeparams) = beststate./c;
  datasummary{f,2} = [beststateinrng,MIrngmin,bestinrngk,...
    totpntltMIthresh,endpntsltMIthresh,iter2conv,nannealpnts,temp2conv];
  infostruct.mu0 = muinit;
  infostruct.MI0 = d(1,MIpasti);
  infostruct.muexp = mubest;
  infostruct.MIexp = bestMI;
  infostruct.T0 = T0;
  infostruct.tau = tempdecay;
  if exist('nappect'); infostruct.naccept = sum(accvec(1:reannealstart(2))); end
  infostruct.mu_num = activeparams;
  infostruct.bounds = [lb./c;ub./c];
  infostruct.scaled_bounds = [lb;ub];
  infostruct.time_hrs = time_hr;
  infostruct.niter = length(accvec);
  infostruct.nannealpnts = nannealpnts;
  datasummary{f,3} = infostruct;
  
  if usemrk
    pltindrng = reannealstart;
  else
    pltindrng = [1,accstatei(end) + 1];
  end
  %=====================================================
  % State Space Visualization
  %=====================================================
  maxMI = max(MIstate);
  minMI = min(MIstate);
  
  MIbinned = ((MIstate - minMI)./(maxMI - minMI)).*ncolbins;
  MIbinned = ceil(MIbinned);
  MIbinned(MIbinned==0) = 1;
  MIvals = ([1/ncolbins:1/ncolbins:1].*(maxMI - minMI) + minMI)';
  MIvals = [minMI;MIvals];
  colors = jet(ncolbins);
  switch plottype
    case '3d'
      %=====================================================
      % PLOT 1: Points plotted in 3-D with color indicating MI
      %=====================================================
      plotlabel = '3d';
      figure;
      for p = 1:ncolbins
        MIindp = find(MIbinned == p);
        MIindabs = accstatei(MIindp);
        for m = 1:length(pltindrng) -1
          MIind = find(MIindabs >= pltindrng(m) & MIindabs < pltindrng(m+1));
          MIind = MIindp(MIind);
          %MIind = find(MIbinned == p);
          if nvar == 2
            plotdata = [accstate(MIind,1),accstate(MIind,2),MIstate(MIind)];
            zlbl = 'Mutual Information';
            zlm = [minMI,maxMI];
          elseif nvar == 3
            plotdata = [accstate(MIind,1),accstate(MIind,2),accstate(MIind,3)];
            zlbl = paramstr{activeparams(3)};
            zlm = [lb(activeparams(3)),ub(activeparams(3))];
          end
          plot3(plotdata(:,1),plotdata(:,2),plotdata(:,3),...
            markers{mod(m,5)+5*0^mod(m,5)},'markerfacecolor',colors(p,:),...
            'markeredgecolor','k',...
            'linewidth',0.5,...
            'displayname',sprintf('%d) %0.2f-%0.2f',m,MIvals(p),MIvals(p+1)));
          hold on;
        end
      end
      hold off; grid on;
      xlabel(paramstr{activeparams(1)}); xlim([lb(activeparams(1)),ub(activeparams(1))]);
      ylabel(paramstr{activeparams(2)}); ylim([lb(activeparams(2)),ub(activeparams(2))]);
      zlabel(zlbl); zlim(zlm);
      title(sprintf('%s',fstr))
      if svp
        set(gcf,'PaperPositionMode','auto');
        saveas(gcf,...
          strcat(path,sprintf('%s_%s_%s.png',fnameheader,plotlabel,fnamefooter)))
      end
    case '2d'
      %=====================================================
      % PLOT 2: 2-D slices with color indicating MI
      %=====================================================
      [a1sort,ind] = sort(accstate(:,1));
      pps = round(npnts/nslices); %pnts per slice
      for n = 1:nslices
        plotlabel = sprintf('2d_slice%02d',n);
        figure; hold on;
        slsp = ind((1+pps*(n-1)):min(pps*n,npnts)); % slice points
        MIbinsls = MIbinned(slsp);
        for p = 1:ncolbins
          MIslsp = find(MIbinsls == p);
          if ~isempty(MIslsp)
            plot(accstate(slsp(MIslsp),2),accstate(slsp(MIslsp),3),...
              'o','markerfacecolor',colors(p,:), 'markeredgecolor','k',...
              'linewidth',0.5,...
              'displayname',sprintf('%0.2f-%0.2f',MIvals(p),MIvals(p+1)));
          end
        end
        title({[fstr],...
          [sprintf('Slice %d, %s = [%0.2f,%0.2f], %d pnt density',n,paramstr{activeparams(1)},...
          a1sort(1+pps*(n-1)),a1sort(min(pps*n,npnts)), length(slsp) )]});
        xlabel(paramstr{activeparams(2)}); xlim([lb(activeparams(2)),ub(activeparams(2))]);
        ylabel(paramstr{activeparams(3)}); ylim([lb(activeparams(3)),ub(activeparams(3))]);
        if svp
          set(gcf,'PaperPositionMode','auto');
          saveas(gcf,...
            strcat(path,sprintf('%s_%s_%s.png',fnameheader,plotlabel,fnamefooter)))
        end
      end % PLOT 2
      
    case 'lin'
      %=====================================================
      % PLOT 3: Plots
      %=====================================================
      figure('Name','Tx Parameters','position',[98   451   783   498]);
      MIindabs = accstatei(:);
      pltindrng = [pltindrng;max(MIindabs)-1];
      if pltindrng(end) < pltindrng(end-1); pltindrng = pltindrng(1:end-1); end
      for m = 1:length(pltindrng) -1
        MIind = find(MIindabs >= pltindrng(m) & MIindabs < pltindrng(m+1));
        plotdata = accstate(MIind,:);
        for v = 1:length(activeparams)
          if m == 1
            hand(v) = subplot(1+(length(activeparams) > 3),3,v);
          end
          
          if isempty(endpntsltMIthresh) || length(endpntsltMIthresh) < m || ...
              endpntsltMIthresh(m) <= 0.5
            reannconv = '';
          elseif endpntsltMIthresh(m) > 0.5
            reannconv = '*';
          end
          
          set(gcf,'currentaxes',hand(v))
          plot(accstate(MIind,v),MIstate(MIind),...
            markers{mod(m,5)+5*0^mod(m,5)},'markerfacecolor',colors(m,:),...
            'markeredgecolor','k',...
            'linewidth',0.5,...
            'displayname',sprintf('ReAnneal %d%s',m,reannconv));
          hold on; grid on;
          axis tight
          %xlim([lb(activeparams(v))./c(activeparams(v)),ub(activeparams(v))./c(activeparams(v))]);
          %ylim([minMI - minMI*0.1,maxMI + minMI*0.1])
          xlabel(paramstr{activeparams(v)});
        end
      end
      for v = 1:length(activeparams)
        set(gcf,'currentaxes',hand(v))
        plot([lb(activeparams(v))./c(activeparams(v)),ub(activeparams(v))./c(activeparams(v))],...
          [maxMI,maxMI],'--r');
        plot([beststate(v)./c(activeparams(v)),beststate(v)./c(activeparams(v))],...
          [minMI - minMI*0.1,maxMI],'--r');
        plot(beststate(v)./c(activeparams(v)),maxMI,'pr');
        hold on;
      end
      %       for p = 1:ncolbins
      %         MIindp = find(MIbinned == p);
      %         MIindabs = accstatei(MIindp);
      %         for m = 1:length(pltindrng) -1
      %           MIind = find(MIindabs >= pltindrng(m) & MIindabs < pltindrng(m+1));
      %           MIind = MIindp(MIind);
      %           %MIind = find(MIbinned == p);
      %           plotdata = accstate(MIind,:);
      %           for v = 1:length(activeparams)
      %             if m == 1
      %               hand(v) = subplot(1,length(activeparams),v);
      %             end
      %             set(gcf,'currentaxes',hand(v))
      %             plot(plotdata(:,v),MIstate(MIind),...  ones(1,length(MIind)), plotdata(:,v),...
      %               markers{mod(m,5)+5*0^mod(m,5)},'markerfacecolor',colors(p,:),...
      %               'markeredgecolor','k',...
      %               'linewidth',0.5,...
      %               'displayname',sprintf('%d) %0.2f-%0.2f',m,MIvals(p),MIvals(p+1)));
      %             hold on; grid on;
      %             xlim([lb(activeparams(v))./c(activeparams(v)),ub(activeparams(v))./c(activeparams(v))]);
      %           end
      %         end
      %       end
      hold off; grid on;
      set(gcf,'currentaxes',hand(2))
      title(sprintf('%s',fstr))
      
      set(gcf,'currentaxes',hand(1)); ylabel('MI');
      set(gcf,'currentaxes',hand(4)); ylabel('MI');
      if svp
        set(gcf,'PaperPositionMode','auto');
        saveas(gcf,...
          strcat(path,sprintf('%s_%s_%s.png',fnameheader,plotlabel,fnamefooter)))
      end
  end % switch plot str
  
  if isempty(Fimg) + isempty(Mimg) == 0
    if isempty(pF) || ismember(f,pF)
      Fedgth = 60; %0.04;%
      if strmatch('Medgeth',PropertyNames)
        Medgth = PropertyVal{strmatch('Medgeth',PropertyNames)};
      else
        Medgth = 20;
      end
      SIth = 100;
      jetmap = jet(256);
      if isstruct(Fimg) || isstruct(Mimg)
        elstr = fname((slashi(end-1)+1):(slashi(end)-1));
        if isstruct(Fimg)
          reF = imresize(getfield(Fimg,elstr),[Fm,Fn]);
          Fth = auto_thresh(reF,'auto',[]);
          reF(reF < Fth(1)) = Fth(1);
          reF(reF > Fth(2)) = Fth(2);
          reM = imresize(Mimg,[Fm,Fn]);
          Fedgth = SIth;
        else
          reM = imresize(getfield(Mimg,elstr),[Fm,Fn]);
          reF = imresize(Fimg,[Fm,Fn]);
          Medgth = SIth;
        end
      else
        reF = imresize(Fimg,[Fm,Fn]);
        reM = Mimg;%imresize(Mimg,[Fm,Fn]);
      end
      reMinit = transform_image(reM,muinit,'extrapval',min(reM(:)),'outputimsz',size(reF)); %imoverlay1x4(F,'',Mp,Mt,'')
      reMbest = transform_image(reM,mubest,'extrapval',min(reM(:)),'outputimsz',size(reF)); %imoverlay1x4(F,'',Mp,Mt,'')
      
      if pMI
        twoimgMIkde(reF,reM,'mu',muinit,'plot',pMI);
        title(sprintf('%s Initial State',fstr))
        twoimgMIkde(reF,reM,'mu',mubest,'plot',pMI);
        title(sprintf('%s Best State',fstr))
      end
      
      imoverlay(reF,'',reMinit,reMinit,'0',varargin{:},'npanels',1)
      title(sprintf('%s Initial State',fstr))
      
      imoverlay(reF,'',reMbest,reMbest,'exp',varargin{:},'npanels',1)
      title(sprintf('%s Best State',fstr))
    end
  end % PLOTS moving and fixed images
  
  dockf on all
  %Mt = transform_image(M,muexp,'extrapval',nan); %imoverlay1x4(F,'',Mp,Mt,'')
  clear infostruct
end % FOR files in Fcell