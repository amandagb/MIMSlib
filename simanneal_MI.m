function [mu,exp_outputs,phiMI] = simanneal_MI(Fin,Min,varargin)
%% Script:
% Description:
% Example:
% Required Functions: twoimgMIkde
% INPUTS ----------------------------------------------------------------
% Fin:  fixed image
% Min:  moving image
% varargin - 'PropertyName','PropertyValue'
%   • 'mu0':    1x6 vector with initial transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
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
%  22 Oct 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  25 Oct 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Changed the simulating annealing function so that only one variable
%     is varied rather than all six
%   8 Nov 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Allowed a variable number of registration parameters
%  15 Nov 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.1
%     Changed parameter factors so that all would be approximately on the
%     order of 10's (scales require 10X factor)
%  16 Jan 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.2
%     Changed parameter factors so all would be on order of 1's (10-20 max)
%  15 Oct 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     Fitting parameters to work with new twoimgMIkde function
%  15 Nov 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     4


global M F muItx nbins klogcell pF varparam mu0 x bestx fval bestfval ...
  temp tau gaussvec gaussvec_unit ...
  c tx_scalefactor ty_scalefactor rot_scalefactor sx_scalefactor sy_scalefactor k_scalefactor
M = Min;
F = Fin;
[Fm,Fn] = size(F);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

%.....................................................................
%   simanneal_MI varargin
%.....................................................................

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

x = zeros(nparam,initvars);
bestx = zeros(nparam,initvars);
fval = zeros(1,initvars);
bestfval = zeros(1,initvars);
temp = zeros(nparam,initvars);

% Set simulated annealing options
optSA = saoptimset('simulannealbnd');%optimoptions(@simulannealbnd);%
optSA.PlotFcns = {@saplotfall,@saplot_temp,@saplotvar,@plotIMoverlay};%,@plotannealupdate};%@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplot_temp,@saplotvar};

if PropertyVal{strmatch('track_gen',PropertyNames)}
  optSA.AnnealingFcn = @track_annealingfast;
end

if PropertyVal{strmatch('track_accept',PropertyNames)}
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
switch MImethod
  case 'mattes'
    
  case 'kde'
    if sumdim
      [mu,phiMI] = simulannealbnd(@computeMIkdesum,mu0(varparam),lb(varparam),ub(varparam),optSA);
    else
      [mu,phiMI] = simulannealbnd(@computeMIkde,mu0(varparam),lb(varparam),ub(varparam),optSA);
    end
end
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
% MI Computation Function
%=====================================================
  function SMI = computeMIkde(x); % would need M, nbins, klog_cell, pF, LF to be global
    mu = mu0;
    mu(varparam) = x;
    mu = mu./c;
    MI = twoimgMIkde(F,M,varargin{:},'mu',mu);%'nbins',nbins);
    %for i = 1:size(F,3)
    %  MI(i) = twoimgMIkde(F(:,:,i),M,varargin{:},'mu',mu);%'nbins',nbins);
    %end
    SMI = -MI;
  end

  function SMI = computeMIkdesum(x); % would need M, nbins, klog_cell, pF, LF to be global
    mu = mu0;
    mu(varparam) = x;
    mu = mu./c;
    if length(size(F)) > 2
      for i = 1:size(F,3)
        if sum(sum(F(:,:,i)))
          MI(i) = twoimgMIkde(F(:,:,i),M,varargin{:},'mu',mu);%'nbins',nbins);
        else
          MI(i) = 0;
        end
      end
    else
      MI = twoimgMIkde(F,M,varargin{:},'mu',mu);%'nbins',nbins);
    end
    SMI = -sum(MI);
  end

  function stop = saplotvar(options, optimvalues,flag)
    stop = false;
    param_names = {'trans x','trans y', 'rot', 'scale x', 'scale y', 'skew'};
    param_abb = {'t_x','t_y', '\theta^\circ', 's_x', 's_y', 'k'};
    switch flag
      case 'init'
        plotvar = plot(optimvalues.iteration,...
          optimvalues.x(min(length(optimvalues.x),varparam)), '.');
        set(plotvar,'Tag','saplotvar');
        hold on;
        
        plotvarbest = plot(optimvalues.iteration,...
          optimvalues.bestx(min(length(optimvalues.x),varparam)), 'x');
        set(plotvarbest,'Tag','saplotvarbest');
        
        xlabel('Iteration','interp','none');
        ylabel(sprintf('Parameter values'),'interp','none')
        title(sprintf('Best point: %g',...
          optimvalues.bestx(min(length(optimvalues.x),varparam))),'interp','none');
      case 'iter'
        titlestr1 = '';
        titlestr2 = '';
        plotvar = findobj(get(gca,'Children'),'Tag','saplotvar');
        for i = 1:length(plotvar)
          newX = [get(plotvar(i),'Xdata') optimvalues.iteration];
          newY = [get(plotvar(i),'Ydata') optimvalues.x(i)];
          set(plotvar(i),'Xdata',newX, 'Ydata',newY);
        end
        hold on;
        
        plotvarbest = findobj(get(gca,'Children'),'Tag','saplotvarbest');
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
        
        %set(get(gca,'Title'),'String',sprintf('Best point: %g',optimvalues.bestx(min(length(optimvalues.x),varparam))));
        set(get(gca,'Title'),'String',{sprintf('$$ %s $$',titlestr1),sprintf('$$ %s $$',titlestr2)},'interpreter','latex');
    end
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
        
        xlabel('Iteration','interp','none');
        ylabel('- MI value','interp','none')
        title(sprintf('Current MI: %1.3g; Best MI: %1.3g',...
          optimvalues.fval,optimvalues.bestfval),'interp','none');
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
        set(get(gca,'Title'),'String',sprintf('Current MI: %1.3g; Best MI: %1.3g',...
          optimvalues.fval,optimvalues.bestfval));
    end
  end

  function stop = saplot_temp(options, optimvalues,flag)
    stop = false;
    param_names = {'trans x','trans y', 'rot', 'scale x', 'scale y', 'skew'};
    switch flag
      case 'init'
        plottemp = plot(optimvalues.iteration,optimvalues.temperature(min(length(optimvalues.x),varparam)), '.b');
        set(plottemp,'Tag','saplot_temp');
        xlabel('Iteration','interp','none');
        ylabel(sprintf('Temperature',param_names{varparam}),'interp','none')
        title(sprintf('Annealing Temperature'),'interp','none');
      case 'iter'
        plottemp = findobj(get(gca,'Children'),'Tag','saplot_temp');
        for i = 1:length(plottemp)
          newX = [get(plottemp(i),'Xdata') optimvalues.iteration];
          newY = [get(plottemp(i),'Ydata') optimvalues.temperature(i)];
          set(plottemp(i),'Xdata',newX, 'Ydata',newY);
        end
        %set(get(gca,'Title'),'String',sprintf('Best Parameter Value: %g',optimvalues.bestx(p)));
    end
  end

  function stop = plotIMoverlay(options, optimvalues,flag)
    stop = false;
    switch flag
      case 'init'
        Tx = muItx; Tx(varparam) = optimvalues.x; Tx = Tx./c;
        Mt = transform_image(M, Tx,'extrapval',nan,'outputimsz',size(F)); %imoverlay1x4(F,'',Mp,Mt,'')
        if length(size(Mt)) > 2
          M1t = sum(Mt.*(1/size(Mt,3)),3);
        else
          M1t = Mt;
        end
        Mtedge = edge(M1t,20);
        reF = (F - min(F(:)))./(max(F(:)) - min(F(:)));
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
        title({[sprintf('[%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f]',optimvalues.x)],...
          [sprintf('MI = %1.2f',optimvalues.fval)]});
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
          Mtedge = edge(M1t,20);
          reF = (F - min(F(:)))./(max(F(:)) - min(F(:)));
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