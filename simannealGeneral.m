function [mu,exp_outputs,phiMI] = simannealGeneral(Fin,Min,varargin)
%% Script: simannealGeneral
% Description:
% Example:
% Required Functions: twoimgMIkde
% INPUTS ----------------------------------------------------------------
% Fin:  M x N x df fixed image (typically taken to be double precision)
% Min:  M x N x dm moving image
% varargin - 'PropertyName','PropertyValue'
%   • 'mu0':    1x6 vector with initial transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%   • 'Fmask'/'Mmask':  M x N logical image indicating region of
%   interest/foreground (1) and background (0)
%   • 'nbins':  Number of histogram bins for joint and marginal histogram
%   approximations
%   • 'paramvar': integer or vector indicate which parameter(s) is being
%   varied
%   • 'bounds': integer indicating the upper bounds for each of the affine
%   parameters. The lower bounds will be derived from these.
%   (DEFAULT = [Fn/2,Fm/2,pi,4,4,pi/4])
%   • 'maxiter': integer specifying the maximum number of iterations allowed
%   [DEFUALT = Inf]
%   • 'tempfactor': integer < 1 which specifies the temperature decay rate. The
%   function will be: T = T0*tempfactor^k where k indicates the iteration.
%   Input 1 for simulanneal's default factor: 0.95 [DEFAULT = 0.985]
%   • 't0': integer indicating the intial temperature [DEFAULT = 3]
%   • 'reannealint': integer indicating the number of iterations for
%   reannealing interval [DEFAULT = Inf]
%   • 'track_gen': binary indicating whether to use track_annealingfast
%   function or not which saves nvar normal random #; unitized random
%   vector; x-increment added to the states
%   • 'track_accept': binary indicating whether or to use track_acceptancesa
%   function or not which saves MIold, MInew, MInew - MIold, h1 (boltz),
%   h2 (exp), Uni(0,1), Temp, accept (0/1), new state
%   • 'costfctn': string or numeric indicating the cost function. IF
%   NUMERIC, then the lp norm is used where the numeric is interpreted as
%   the value of p
%   ** NOTE: one can indicate computation of a normalized version of the
%   cost function where the cost is the average value of the included
%   pixels (that is 1/(# non-nan pixels included))
%     >>> For normalized version, prefix string with an "n" <<< ('nmi')
%       OPTIONS:  'mi' (mutual information) [DEFAULT]
%                 'cc' (cross correlation)
%                 'div' (Csizar's I-divergence)
%                 'tv' (total variation),
%                 'max1' (Genab and Reynolds)
%                 'HL' (Herbert/Lechy)
%                 'Huber'
%                 'truncL2' (truncated L2)
%   • 'normcostf': logical indicating whether cost function will be
%   normalized by # of pixels in "true" region of mask (1) or not (0)
%   [DEFAULT = 0]
%   • 'plotprogress': binary indicating whether to plot the output
%   real-time (1) or not (0) [DEFAULT = 1]
%   ..................................................................
%   twoimgMIkde varargin
%   ..................................................................
%   • 'mu': (DEFAULT = [0,0,0,1,1,0] = [tx,ty,rad,sx,sy,sk])
%   • 'nbins': [DEFAULT = 50]
%   • 'logbase': [DEFUALT = 'e']
%   • 'pA'   • 'Abinlim'    • 'axh'   • 'Bbinlim'
%   • 'plot': [DEFAULT = 0]
%   • 'kernel': {'string',bandwidth,dimensions} [DEFAULT = {'norm',[]}]
%             {'kernel1',[bw_d1,bw_d2,...],[d1,d2,...];
%             'kernel2',[bw_d1,bw_d2,...],[d1,d2,...]}
%   ..................................................................
%   transform_image varargin
%   ..................................................................
%   • 'extrapval': [DEFAULT = nan]
%   • 'interpmethod': [DEFAULT = 'bicubic']
%   • 'transtype': [DEFAULT = 'physical']
%   • 'fixedstr'  • 'fixedvals'   • 'varstr'    • 'varvals'
%   • 'rotpnt': [DEFAULT = 'center']
%   • 'shape': [DEFAULT = 'same']
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%   8 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1
%     Based on v4 of simanneal_MI
%  13 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1.1
%     Added scaling of cost function
%  21 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1.2
%     call computeCosts rather than switching
%  28 May 2015   AGBalderrama      amandagbalderrama@gmail.com     1.3
%     Automated scaling for cost functions based on initial cost to be
%     equal to exactly 200

global M Mm Mn Md F Fvec Fm Fn Fd muItx Fmask Mmask Fnan Mnan...
  normcost nbins varparam mu0 x bestx fval bestfval ...
  temp tau gaussvec gaussvec_unit pval...
  c ub lb scaleJ costfstr fullcostname maxiter tx_scalefactor ty_scalefactor rot_scalefactor sx_scalefactor sy_scalefactor k_scalefactor
M = Min;
[Mm,Mn,Md] = size(M);

F = Fin;
[Fm,Fn,Fd] = size(F);
Fvec = reshape(F(:),[Fm*Fn,Fd]);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

%.....................................................................
%   simanneal_MI varargin
%.....................................................................

if strmatch('Fmask',PropertyNames)
  Fmask = PropertyVal{strmatch('Fmask',PropertyNames)};
else
  Fmask = true(Fm,Fn);
end

if strmatch('Mmask',PropertyNames)
  Mmask = PropertyVal{strmatch('Mmask',PropertyNames)};
else
  Mmask = true(Mm,Mn);
end

if strmatch('paramvar',PropertyNames)
  varparam = PropertyVal{strmatch('paramvar',PropertyNames)};
else
  varparam = [1:6]; %tx, ty, theta, sx, sy, sk
end
nparam = length(varparam);

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 50;
end

Fnan = F;
Fnan(~repmat(Fmask,[1,1,Fd])) = nan;

Mnan = M;
Mnan(~repmat(Mmask,[1,1,Md])) = nan;

%=====================================================
% Simulated Annealing Parameters
%=====================================================
c = sastatescalefactors;
muItx = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk

if strmatch('mu0',PropertyNames)
  mu0 = PropertyVal{strmatch('mu0',PropertyNames)};
else
  mu0 = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk
end
mu0 = mu0.*c;
def_bnd = [Fn*0.1,Fm*0.1,10*pi/180,1.3,1.3,0.2].*c;
ub = repmat(muItx,2,1);
if strmatch('bounds',PropertyNames)
  ub = PropertyVal{strmatch('bounds',PropertyNames)};
  ub = ub(:,varparam);
  ub = ub.*repmat(c(:,varparam),size(ub,1),1);
else
  ub = mu0 + def_bnd;
  ub(4:6) = def_bnd(4:6);
end
if size(ub,1) > 1
  lb = ub(2,:);
  ub = ub(1,:);
else
  lb = mu0 - def_bnd;
  lb(4:5) = 1./def_bnd(4:5).*c(4:5).^2;
  lb(6) = -def_bnd(6);
end

if strmatch('maxiter',PropertyNames)
  maxiter = PropertyVal{strmatch('maxiter',PropertyNames)};
  initvars = maxiter;
else
  maxiter = Inf;
  initvars = 3000;
end

if strmatch('tempfact',PropertyNames)
  tau = PropertyVal{strmatch('tempfact',PropertyNames)};
else
  tau = 0.985;
end

if strmatch('t0',PropertyNames)
  t0 = PropertyVal{strmatch('t0',PropertyNames)};
else
  t0 = 3;
end

if strmatch('plot',PropertyNames)
  viewprogress = PropertyVal{strmatch('plot',PropertyNames)};
else
  viewprogress = 1;
end

x = zeros(nparam,initvars);
bestx = zeros(nparam,initvars);
fval = zeros(1,initvars);
bestfval = zeros(1,initvars);
temp = zeros(nparam,initvars);

% Set simulated annealing options
optSA = saoptimset('simulannealbnd');
if viewprogress
  optSA.PlotFcns = {@saplotfall,@saplot_temp,@saplotvar,@overlayCH};%@overlayGrayEdge};%,@plotannealupdate};%@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplot_temp,@saplotvar};
else
  optSA.PlotFcns = {};
end

if ~isempty(PropertyVal(strmatch('track_gen',PropertyNames)))
  optSA.AnnealingFcn = @track_annealingfast;
end

if ~isempty(PropertyVal(strmatch('track_accept',PropertyNames)))
  optSA.AcceptanceFcn = @track_acceptancesa;
end

if strmatch('reann',PropertyNames)
  reanint = PropertyVal{strmatch('reann',PropertyNames)};
  if strmatch(reanint,'scale')
    reanint = round((log(1e-2/t0)/log(tau))/2*max(1,nparam/4));
  end
else
  reanint = Inf;
end

if strmatch('MImethod',PropertyNames)
  MImethod = PropertyVal{strmatch('MImethod',PropertyNames)};
else
  MImethod = 'kde';
end

if strmatch('costf',PropertyNames)
  costfstr = PropertyVal{strmatch('costf',PropertyNames)};
else
  costfstr = 'mi';
end

if strmatch('pval',PropertyNames)
  pval = PropertyVal{strmatch('pval',PropertyNames)};
else
  pval = 1; %tx, ty, theta, sx, sy, sk
end

if strmatch('normc',PropertyNames)
  normcost = PropertyVal{strmatch('normc',PropertyNames)};
else
  normcost = 0;
end

if strmatch('sumdim',PropertyNames)
  sumdim = PropertyVal{strmatch('sumdim',PropertyNames)};
else
  sumdim = 0;
end

if ~isempty(tau)% ~= 1
  optSA.TemperatureFcn = @customtemp;
end
optSA.InitialTemperature = t0;
optSA.ReannealInterval = Inf;
optSA.StallIterLimit = Inf; %1000*nparam;
optSA.MaxFunEvals = Inf;
optSA.OutputFcns = @collectvars;
optSA.MaxIter = maxiter;
% optSA.TimeLimit = 1800;
optSA.ReannealInterval = reanint;
switch costfstr
  case 'mi'
    Jrng = [-1e-1,0];
    fullcoststr = '- MI';
    normcost = 0;
  case 'mine'
    Jrng = [-log(max(sum(Fmask(:)),sum(Mmask(:)))),0];
    fullcoststr = 'MIC';
    normcost = 0;
  case 'lp'
    Jrng = [0,abs(sum((nanmax(abs(Fnan(:))) - nanmin(abs(Fnan(:)))).^pval))];
    fullcoststr = sprintf('$$\\ell%1.2f$$',pval);
  case 'lcca'
    Jrng = [0,5];%[0,max(sum(Fmask(:)),sum(Mmask(:)))];
    fullcoststr = sprintf('Hein LCCA');
  case 'cca'
    Jrng = [0,1e-5];%[0,min(Fd,Md)];
    fullcoststr = sprintf('MATLAB CCA');
  case 'div'
    Jrng = [0,log(max(sum(Fmask(:)),sum(Mmask(:))))];
    fullcoststr = sprintf('$$D_{KL}(p\\|n^{-1})$$');
  case 'tv'
    Jrng = [0,30];%[abs(sum((nanmean(abs(Fnan(:))) - nanmean(abs(Mnan(:)))).^pval)),...
      %abs(sum((nanmax(abs(Fnan(:))) - nanmean(abs(Mnan(:)))).^pval))];
    fullcoststr = sprintf('TV');
    normcost = 0;
  case 'max1'
    Jrng = [0,1];
    fullcoststr = 'max1';%sprintf('$$d_{\\ell%1.2f}(d_{\\ell%1.2f} + 1)^{-1}$$',pval,pval);
  case 'HL'
    Jrng = log(sum((nanmax(abs(Fnan(:))) - nanmin(abs(Fnan(:)))).^pval).*[-1,1]);
    fullcoststr = sprintf('$$log(d_{\\ell%1.2f} + 1)$$',pval);
end
[J,totpxls] = computeCost(F,M,costfstr,'mu0',mu0./c,'Fmask',Fmask,'Mmask',Mmask,'pval',pval);
targetJ0 = 200;
scaleJ = log10((targetJ0*totpxls^(0^abs(normcost-1)))/J);% (floor(log10(t0))+1) - (log10(abs(diff(Jrng)))-1) - log10(sum(Fmask(:)))*0^normcost + 1;
fullcostname = sprintf('%s $$\\cdot 10^{%1.2f}$$',fullcoststr,scaleJ);
[mu,phiMI] = simulannealbnd(@computeCostsub,mu0(varparam),lb(varparam),ub(varparam),optSA);

muout = [0,0,0,1,1,0].*c;
muout(varparam) = mu;
muout = muout./c;
% muout(3) = muout(3)/rot_scalefactor;
% muout([4,5]) = muout([4,5])/sx_scalefactor;
% muout(6) = muout(6)/k_scalefactor;
mu = muout;

exp_outputs.x = x;
exp_outputs.bestx = bestx;
exp_outputs.fval = fval;
exp_outputs.bestfval = bestfval;
exp_outputs.temp = temp;
exp_outputs.nbins = nbins;
%=====================================================
% Cost Function Computation
%=====================================================
  function J = computeCostsub(x)%,costfstr,varargin)
    mu = mu0;
    mu(varparam) = x;
    mu = mu./c;
    TxM = transform_image(M,mu,'extrapval',M(1));
    TxMask = transform_image(double(Mmask),mu,'extrapval',double(Mmask(1)));
    usei = and(Fmask,TxMask);
    TxMvec = reshape(TxM(:),[Mm*Mn,Md]);
    if Fd > Md
      dFM = Fvec - [TxMvec,zeros(Fm*Fn,Fd - Md)];
    elseif Md > Fd
      dFM = [Fvec,zeros(Fm*Fn,Md - Fd)] - TxMvec;
    else
      dFM = Fvec - TxMvec;
    end
    
    switch costfstr
      case 'mi'
        switch MImethod
          case 'mattes'
            
          case 'kde'
            %Jrng = [
            %if sumdim
            % if Fd > 1
            %    for i = 1:Fd
            %      if sum(sum(Fnan(:,:,i)))
            %        MI(i) = twoimgMIkde(Fnan(:,:,i),Mnan,varargin{:},'mu',mu);%'nbins',nbins);
            %      else
            %        MI(i) = 0;
            %      end
            %    end
            %  else
            %    MI = twoimgMIkde(Fnan,Mnan,varargin{:},'mu',mu);%'nbins',nbins);
            %  end
            %  J = -sum(MI);
            %else
            [MI,outinfo] = twoimgMIkde(Fnan,Mnan,varargin{:},'mu',mu);%'nbins',nbins);
            totpxls = outinfo.nonNANpxls;
            J = -MI;
            %end
        end
      case 'lp'
        %dvec = Fvec(usei(:),:) - TxMvec(usei(:),:);
        %J = norm(dvec(:),pval); --- use of matlab native function could be more efficient
        dFM = sum(abs(dFM).^pval,2);
        Jimg = reshape(dFM,[Fm,Fn]);
        J = sum(Jimg(usei));
        if normcost
          J = J/sum(usei(:));
        end
      case 'mine' % Albanese, minepy-1.0.0, Bioinformatics (2013) see: minepy.sourceforge.net (based on Reshef 2011 MIC metric paper)
        mineout = mine(Fvec(usei(:),:)',TxMvec(usei(:),:)');
        J = mineout.mic;
      case 'lcca' % Heinrich, MICCAI 2014, see http://www.mpheinrich.de/software.html
        Jimg = LCCA_2d(F, TxM, 5, 1e-3);
        J = sum(Jimg(usei));
        if normcost
          J = J/sum(usei(:));
        end
      case 'cca'
        %[A,B,r,U,V,stats] = canoncorr(Fvec(usei(:),:),TxMvec(usei(:),:));
        [~,~,J,~,~,~] = canoncorr(Fvec(usei(:),:),TxMvec(usei(:),:));
        J = sum(J);
        if normcost
          J = J/sum(usei(:));
        end
      case 'div' % J(x,x0) = sum_i x_i log(x_i/x_0i) + sum_i (x_i - x_0i), normalize: sum_i x_i = 1, uniform prior: x_0i = 1/N => J(x) = sum_i x_i log x_i
        %x = abs(dFM(usei))./sum(abs(dFM(usei)));
        J = entropy(dFM(usei))*log(2); % in NATs
      case 'tv'
        dFM = sum(abs(dFM).^pval,2);
        Jimg = reshape(dFM,[Fm,Fn]);
        J = nanmax(dFM) - nanmin(dFM);
      case 'max1'
        dFM = sum(abs(dFM).^pval,2);
        Jimg = reshape(dFM./(1+dFM),[Fm,Fn]);
        J = sum(Jimg(usei));
        if normcost
          J = J/sum(usei(:));
        end
      case 'HL'
        dFM = sum(log(abs(dFM).^pval + 1),2);
        Jimg = reshape(dFM,[Fm,Fn]);
        J = sum(Jimg(usei));
        if normcost
          J = J/sum(usei(:));
        end
    end
    
    J = J*10^(scaleJ);
  end

  function stop = saplotvar(options, optimvalues,flag)
    stop = false;
    param_names = {'trans x','trans y', 'rot', 'scale x', 'scale y', 'skew'};
    param_abb = {'t_x','t_y', '\theta', 's_x', 's_y', 'k'};
    titlestr1 = '';
    titlestr2 = '';
    linecolrs = lines(6);
    switch flag
      case 'init'
        plotvar = plot(optimvalues.iteration,...
          optimvalues.x(min(length(optimvalues.x),varparam)), '.');
        set(plotvar,'Tag','saplotvar');
        hold on;
        
        plotvarbest = plot(optimvalues.iteration,...
          optimvalues.bestx(min(length(optimvalues.x),varparam)), 'x');
        set(plotvarbest,'Tag','saplotvarbest');
        
        xlabel('Iteration ($$i$$)','interp','latex');
        ylabel(sprintf('Parameter values'),'interp','none')
        for i = 1:length(plotvarbest)
          bestX = [get(plotvarbest(i),'Xdata') optimvalues.iteration];
          bestY = [get(plotvarbest(i),'Ydata') optimvalues.bestx(i)];
          set(plotvarbest(i),'Xdata',bestX, 'Ydata',bestY);
          
          if i > 3
            titlestr2 = sprintf('%s %s = %1.4f;',titlestr2,param_abb{varparam(i)},optimvalues.bestx(i)./c(varparam(i)));
          else
            titlestr1 = sprintf('%s %s = %1.4f;',titlestr1,param_abb{varparam(i)},optimvalues.bestx(i)./c(varparam(i)));
          end
        end
        title({sprintf('$$ %s $$',titlestr1),sprintf('$$ %s $$',titlestr2)},'interpreter','latex');
      case 'iter'
        plotvar = sort(findobj(get(gca,'Children'),'Tag','saplotvar'));
        for i = 1:length(plotvar)
          newX = [get(plotvar(i),'Xdata') optimvalues.iteration];
          newY = [get(plotvar(i),'Ydata') optimvalues.x(i)];
          set(plotvar(i),'Xdata',newX, 'Ydata',newY);
        end
        hold on;
        
        plotvarbest = sort(findobj(get(gca,'Children'),'Tag','saplotvarbest'));
        for i = 1:length(plotvarbest)
          bestX = [get(plotvarbest(i),'Xdata') optimvalues.iteration];
          bestY = [get(plotvarbest(i),'Ydata') optimvalues.bestx(i)];
          set(plotvarbest(i),'Xdata',bestX, 'Ydata',bestY);
          
          if i > 3
            titlestr2 = sprintf('%s %s = %1.4f;',titlestr2,param_abb{varparam(i)},optimvalues.bestx(i)./c(varparam(i)));
          else
            titlestr1 = sprintf('%s %s = %1.4f;',titlestr1,param_abb{varparam(i)},optimvalues.bestx(i)./c(varparam(i)));
          end
        end
        set(get(gca,'Title'),'String',{sprintf('$$ %s $$',titlestr1),sprintf('$$ %s $$',titlestr2)},'interpreter','latex');
        
        if (maxiter - optimvalues.iteration) == 1
          for i = 1:length(plotvar)
            text(optimvalues.iteration+10*i,optimvalues.bestx(i),...
              param_abb{varparam(i)},'color',linecolrs(i,:));
          end
        end
    end
    %     if ~mod(optimvalues.iteration,499) && optimvalues.iteration
    %       blki = (optimvalues.iteration/499);
    %       plot(repmat([blki-1,blki].*499,[length(ub),1])',...
    %         repmat(lb',[1,2])','--')
    %       plot(repmat([blki-1,blki].*499,[length(ub),1])',...
    %         repmat(ub',[1,2])','--')
    %
    %       %       switch rem(blki,3)
    %       %         case 1
    %       %           ind = 1:floor(length(lb)/3);
    %       %         case 2
    %       %           ind = ceil(length(lb)/3) + [1:floor(length(lb)/3)];
    %       %         case 0
    %       %           ind = (ceil(length(lb)/3)+floor(length(lb)/3)):length(lb);
    %       %       end
    %       %       plot(repmat([blki-1,blki].*199,[length(ind),1])',...
    %       %         repmat(lb(ind)',[1,2])','--')
    %       %       plot(repmat([blki-1,blki].*199,[length(ind),1])',...
    %       %         repmat(ub(ind)',[1,2])','--')
    %     end
  end

  function stop = saplotfall(options,optimvalues,flag)
    stop = false;
    %plot the current function value against the iteration number
    switch flag
      case 'init'
        plotf = plot(optimvalues.iteration,optimvalues.fval, '.b');
        set(plotf,'Tag','saplotf');
        hold on;
        
        plotBest = plot(optimvalues.iteration,optimvalues.bestfval, 'xr');
        set(plotBest,'Tag','saplotfbest');
        
        xlabel('Iteration ($$i$$)','interp','latex');
        ylabel(sprintf('J = %s',fullcostname),'interp','latex')%ylabel('- MI value','interp','none')
        title(sprintf('$$J_{%d}$$: %1.3g; $$J_{i^*}$$: %1.3g',...
          optimvalues.iteration,optimvalues.fval,optimvalues.bestfval),'interp','latex');
      case 'iter'
        plotf = findobj(get(gca,'Children'),'Tag','saplotf');
        newX = [get(plotf,'Xdata') optimvalues.iteration];
        newY = [get(plotf,'Ydata') optimvalues.fval];
        set(plotf,'Xdata',newX, 'Ydata',newY);
        hold on;
        plotBest = findobj(get(gca,'Children'),'Tag','saplotfbest');
        bestX = [get(plotBest,'Xdata') optimvalues.iteration];
        bestY = [get(plotBest,'Ydata') optimvalues.bestfval];
        set(plotBest,'Xdata',bestX, 'Ydata',bestY);
        set(get(gca,'Title'),'String',sprintf('$$J_{%d}$$: %1.3g; $$J_{i^*}$$: %1.3g',...
          optimvalues.iteration,optimvalues.fval,optimvalues.bestfval));
    end
    grid on
  end

  function stop = saplot_temp(options, optimvalues,flag)
    stop = false;
    param_names = {'trans x','trans y', 'rot', 'scale x', 'scale y', 'skew'};
    switch flag
      case 'init'
        plottemp = plot(optimvalues.iteration,optimvalues.temperature(min(length(optimvalues.x),varparam)), '.b');
        set(plottemp,'Tag','saplot_temp');
        xlabel('Iteration ($$i$$)','interp','latex');
        ylabel(sprintf('Temperature',param_names{varparam}),'interp','none')
        title(sprintf('Annealing Temperature (\\tau = %1.2f)',tau));
      case 'iter'
        plottemp = findobj(get(gca,'Children'),'Tag','saplot_temp');
        for i = 1:length(plottemp)
          newX = [get(plottemp(i),'Xdata') optimvalues.iteration];
          newY = [get(plottemp(i),'Ydata') optimvalues.temperature(i)];
          set(plottemp(i),'Xdata',newX, 'Ydata',newY);
        end
        %set(get(gca,'Title'),'String',sprintf('Best Parameter Value: %g',optimvalues.bestx(p)));
    end
    grid on
  end

  function stop = overlayGrayEdge(options, optimvalues,flag)
    stop = false;
    switch flag
      case 'init'
        Tx = muItx; Tx(varparam) = optimvalues.x; Tx = Tx./c;
        Mt = transform_image(M, Tx,'extrapval',nan,'outputimsz',[Fm,Fn]); %imoverlay1x4(F,'',Mp,Mt,'')
        if length(size(Mt)) > 2
          M1t = sum(Mt.*(1/size(Mt,3)),3);
        else
          M1t = Mt;
        end
        %>>>> Specific to mean variance equalize images
        M1t = Mt(:,:,1);
        M1t(M1t < -3) = -3;
        M1t(M1t > 3) = 3;
        %<<<<
        Mtedge = edge(M1t,max(M1t(:))*0.10);%20);
        
        %>>>> Specific to mean variance equalize images
        reF = F; reF(reF < -3) = -3; reF(reF > 3) = 3;
        %<<<<
        reF = (reF - min(reF(:)))./(max(reF(:)) - min(reF(:)));
        if length(size(F)) > 2
          ch1 = reF(:,:,1); ch2 = reF(:,:,2);
          if size(F,3) >= 3
            ch3 = reF(:,:,3);
          else
            ch3 = zeros(size(F,1),size(F,2));
          end
          v1 = 1; v2 = 1; v3 = 1;
        else
          ch1 = reF; ch2 = reF; ch3 = reF;
          v1 = 1; v2 = 0; v3 = 0;
        end
        ch1(find(Mtedge == 1)) = v1;
        ch2(find(Mtedge == 1)) = v2;
        ch3(find(Mtedge == 1)) = v3;
        Fimg3 = repmat(ch1,[1,1,3]);
        Fimg3(:,:,1) = ch1;
        Fimg3(:,:,2) = ch2;
        Fimg3(:,:,3) = ch3;
        Itx = image(Fimg3);
        set(Itx,'Tag','saItx')
        axis image; axis ij;
        title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.x./c)],...
          [sprintf('MI = %1.2f',optimvalues.fval)]});
      case 'iter'
        if mod(optimvalues.iteration,50) == 0 || (maxiter - optimvalues.iteration) <= 2
          Itx = findobj(get(gca,'Children'),'Tag','saItx');
          Tx = muItx; Tx(varparam) = optimvalues.x; Tx = Tx./c;
          Mt = transform_image(M, Tx,'extrapval',nan,'outputimsz',size(F)); %imoverlay1x4(F,'',Mp,Mt,'')
          if length(size(Mt)) > 2
            M1t = sum(Mt.*(1/size(Mt,3)),3);
          else
            M1t = Mt;
          end
          %>>>> Specific to mean variance equalize images
          M1t = Mt(:,:,1);
          M1t(M1t < -3) = -3;
          M1t(M1t > 3) = 3;
          %<<<<
          Mtedge = edge(M1t,max(M1t(:))*0.10);%20);
          
          %>>>> Specific to mean variance equalize images
          reF = F; reF(reF < -3) = -3; reF(reF > 3) = 3;
          %<<<<
          reF = (reF - min(reF(:)))./(max(reF(:)) - min(reF(:)));
          if length(size(F)) > 2
            ch1 = reF(:,:,1); ch2 = reF(:,:,2);
            if size(F,3) >= 3
              ch3 = reF(:,:,3);
            else
              ch3 = zeros(size(F,1),size(F,2));
            end
            v1 = 1; v2 = 1; v3 = 1;
          else
            ch1 = reF; ch2 = reF; ch3 = reF;
            v1 = 1; v2 = 0; v3 = 0;
          end
          ch1(find(Mtedge == 1)) = v1;
          ch2(find(Mtedge == 1)) = v2;
          ch3(find(Mtedge == 1)) = v3;
          Fimg3 = repmat(ch1,[1,1,3]);
          Fimg3(:,:,1) = ch1;
          Fimg3(:,:,2) = ch2;
          Fimg3(:,:,3) = ch3;
          image(Fimg3);
          axis image; axis ij;
          title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.x./c)],...
            [sprintf('MI = %1.2f',optimvalues.fval)]});
        end
    end
  end

  function stop = overlayCH(options, optimvalues,flag)
    stop = false;
    switch flag
      case 'init'
        Tx = muItx; Tx(varparam) = optimvalues.x; Tx = Tx./c;
        Mt = transform_image(M, Tx,'extrapval',M(1),'outputimsz',[Fm,Fn]); %imoverlay1x4(F,'',Mp,Mt,'')
        if length(size(Mt)) > 2
          M1t = sum(Mt.*(1/size(Mt,3)),3);
        else
          M1t = Mt;
        end
        %>>>> Specific to mean variance equalize images
        M1t = Mt(:,:,1); M1t(M1t < -3) = -3; M1t(M1t > 3) = 3;
        M1t = (M1t + 3)./6;
        %<<<<
        
        %>>>> Specific to mean variance equalize images
        reF = F; reF(reF < -3) = -3; reF(reF > 3) = 3;
        %<<<<
        reF = (reF - min(reF(:)))./(max(reF(:)) - min(reF(:)));
        Fimg3 = zeros(Fm,Fn,3);
        Fimg3(:,:,1) = M1t;
        %Fimg3(:,:,2) = reF(:,:,1);
        Fimg3(:,:,3) = reF(:,:,1);
        Itx = image(Fimg3);
        set(Itx,'Tag','saItx')
        axis ij;axis tight;%axis image; 
        title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.x./c)],...
          [sprintf('$$J_{%d}$$ = %1.2f',optimvalues.iteration, optimvalues.fval)]},...
          'interpreter','latex');
      case 'iter'
        if mod(optimvalues.iteration,50) == 0
          Itx = findobj(get(gca,'Children'),'Tag','saItx');
          Tx = muItx; Tx(varparam) = optimvalues.x; Tx = Tx./c;
          Mt = transform_image(M, Tx,'extrapval',nan,'outputimsz',size(F)); %imoverlay1x4(F,'',Mp,Mt,'')
          if length(size(Mt)) > 2
            M1t = sum(Mt.*(1/size(Mt,3)),3);
          else
            M1t = Mt;
          end
          %>>>> Specific to mean variance equalize images
          M1t = Mt(:,:,1); M1t(M1t < -3) = -3; M1t(M1t > 3) = 3;
          M1t = (M1t + 3)./6;
          %<<<<
          
          %>>>> Specific to mean variance equalize images
          reF = F; reF(reF < -3) = -3; reF(reF > 3) = 3;
          %<<<<
          reF = (reF - min(reF(:)))./(max(reF(:)) - min(reF(:)));
          Fimg3 = zeros(Fm,Fn,3);
          Fimg3(:,:,1) = M1t;
          %Fimg3(:,:,2) = reF(:,:,1);
          Fimg3(:,:,3) = reF(:,:,1);
          Itx = image(Fimg3);
          set(Itx,'Tag','saItx')
          axis ij;axis tight;%axis image; 
          title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.x./c)],...
            [sprintf('$$J_{%d}$$ = %1.2f',optimvalues.iteration, optimvalues.fval)]},...
            'interpreter','latex');
        elseif (maxiter - optimvalues.iteration) <= 1
          Itx = findobj(get(gca,'Children'),'Tag','saItx');
          Tx = muItx; Tx(varparam) = optimvalues.bestx; Tx = Tx./c;
          Mt = transform_image(M, Tx,'extrapval',nan,'outputimsz',size(F)); %imoverlay1x4(F,'',Mp,Mt,'')
          if length(size(Mt)) > 2
            M1t = sum(Mt.*(1/size(Mt,3)),3);
          else
            M1t = Mt;
          end
          %>>>> Specific to mean variance equalize images
          M1t = Mt(:,:,1); M1t(M1t < -3) = -3; M1t(M1t > 3) = 3;
          M1t = (M1t + 3)./6;
          %<<<<
          
          %>>>> Specific to mean variance equalize images
          reF = F; reF(reF < -3) = -3; reF(reF > 3) = 3;
          %<<<<
          reF = (reF - min(reF(:)))./(max(reF(:)) - min(reF(:)));
          Fimg3 = zeros(Fm,Fn,3);
          Fimg3(:,:,1) = M1t;
          %Fimg3(:,:,2) = reF(:,:,1);
          Fimg3(:,:,3) = reF(:,:,1);
          Itx = image(Fimg3);
          set(Itx,'Tag','saItx')
          axis ij;axis tight;%axis image; 
          title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.bestx./c)],...
            [sprintf('$$J_{i^*}$$ = %1.2f', optimvalues.bestfval)]},...
            'interpreter','latex');
        end
    end
  end

  function stop = plotannealupdate(options, optimvalues,flag)
    stop = false;
    switch flag
      case 'init'
        nel = numel(optimvalues.x);
        saupdate = plot(optimvalues.iteration,zeros(1,nel), '.');
        set(saupdate,'Tag','saupdate');
        hold on;
        
        saupdateN = plot(optimvalues.iteration,zeros(1,nel), 'x');
        set(saupdateN,'Tag','saupdateN');
        
        xlabel('Iteration','interp','none');
        ylabel(sprintf('Var Amount'),'interp','none')
        title(sprintf('$$ x_{t+1} = x_{t} + N(0,1) T_t/|y| $$'),...
          'interpreter','latex');
      case 'iter'
        saupdate = findobj(get(gca,'Children'),'Tag','saupdate');
        for i = 1:length(saupdate)
          newX = [get(saupdate(i),'Xdata') optimvalues.iteration];
          newY = [get(saupdate(i),'Ydata') gaussvec(i)*optimvalues.temperature(i)];
          set(saupdate(i),'Xdata',newX, 'Ydata',newY);
        end
        hold on;
        
        saupdateN = findobj(get(gca,'Children'),'Tag','saplotvarbest');
        for i = 1:length(saupdateN)
          bestX = [get(saupdateN(i),'Xdata') optimvalues.iteration];
          bestY = [get(saupdateN(i),'Ydata') gaussvec_unit(i)*optimvalues.temperature(i)];
          set(saupdateN(i),'Xdata',bestX, 'Ydata',bestY);
        end
    end
  end

  function temperature = customtemp(optimvalues,options)
    temperature = options.InitialTemperature.*tau.^optimvalues.k;
  end

  function [stop,options,optchanged] = collectvars(options,optimvalues,flag)
    stop = false;
    optchanged = false;
    indx = optimvalues.iteration+1;
    x(:,indx) = optimvalues.x;
    bestx(:,indx) = optimvalues.bestx;
    fval(indx) = optimvalues.fval;
    bestfval(indx) = optimvalues.bestfval;
    temp(:,indx) = optimvalues.temperature;
  end

end