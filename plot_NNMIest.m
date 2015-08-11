function outputs = plot_NNMIest(A,m,totsamp,varargin)
%% Script:
% Description:
% Example:  plot_coarsereg('MI',phi_MI,param_cell,param_names,'normalize',1);
% Required Functions:
% INPUTS ----------------------------------------------------------------
% A : MC x S x 4 matrix where the 1st matrix = H(X); 2nd matrix = H(Y); 3rd
% matrix = H(X,Y); and 4th = MI. Data generated in MADLABmtng_20130701 &
% MADLABmtng_20130801 scripts.
% m : vector of samples used for generating the pdf estimation
% totsamp : total number of samples used
% varargin - 'PropertyName','PropertyValue'
%   {1} theorout : 1 x 4 vector with theoretical outputs [DEFAULT = []]
%   • 'metric' : integer (or vector) indicating the metric the user would
%   like to plot. 1 = H(X), 2 = H(Y), 3 = H(X,Y), 4 = MI. [DEFAULT = 1:4]
%   • 'as%' : binary indicating whether to plat as a percent or not [DEFAULT = 1]
%   • 'title': string indicating plot title
%   • 'theorfrac' : fraction of theoretical value to limit y [DEFAULT = 2*max_standard_dev]
%   • 'subplot' : 0 for a figure per metric, 1 for a 1 x 4 subplot, 2 for a 2
%   x 2 subplot
%   • 'ylim': rv x 2 row vector with [MIN y, MAX y] values for each of the
%   four axes
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%   2 Aug  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  21 Jan  2014   Amanda Gaudreau   amanda.gaudreau@gmail.com     2

theorout = varargin{1};
if isstr(theorout)
  theorout = [];
  PropertyNames = varargin(1:2:length(varargin));
  PropertyVal = varargin(2:2:length(varargin));
else
  PropertyNames = varargin(2:2:length(varargin));
  PropertyVal = varargin(3:2:length(varargin));
end

if strmatch('metric',PropertyNames)
  RV = PropertyVal{strmatch('metric',PropertyNames)}; % plot percent
else
  RV = 1:4;
end
RVstr = {'H(X)','H(Y)','H(X,Y)','MI(X,Y)'};

if strmatch('as%',PropertyNames)
  pp = PropertyVal{strmatch('as%',PropertyNames)}; % plot percent
else
  pp = 1;
end

if strmatch('ylim',PropertyNames)
  ya = PropertyVal{strmatch('ylim',PropertyNames)}; % plot percent
else
  ya = [];
end

if pp
  x = m./totsamp.*100;
  xlstr = sprintf('%% sample points (m, %1.2e total)',totsamp);
else
  x = m;
  xlstr = sprintf('sample points');
end

if strmatch('title',PropertyNames)
  titlestr = PropertyVal{strmatch('title',PropertyNames)}; % plot percent
else
  titlestr = '';
end

if strmatch('theorfrac',PropertyNames)
  frac = PropertyVal{strmatch('theorfrac',PropertyNames)}; % plot percent
else
  frac = [0.003,0.003,0.03,0.2];
end

if strmatch('subplot',PropertyNames)
  nsubplotrows = PropertyVal{strmatch('subplot',PropertyNames)}; % plot percent
else
  nsubplotrows = 0;
end

maxMC = size(A,1);
S = size(A,2);
i = 0;
for rv = RV
  i = i+1;
  if ~ nsubplotrows
    figure; hold on;
  elseif i == 1
    figure; 
    subplot(nsubplotrows,4/nsubplotrows,rv);
    hold on;
  else
    subplot(nsubplotrows,4/nsubplotrows,rv); hold on;
  end
  P = get(gca,'position');
  set(gca,'position',[P(1) 0.1475 P(3) 0.7125]);
  plot(x,A(:,:,rv),'b.');
  mustd = [mean(A(:,:,rv));std(A(:,:,rv))];
  mu = mustd(1,:);
  pmstd = [mustd(1,:) + mustd(2,:); mustd(1,:) - mustd(2,:)];
  plot(x,mu,'.r');%,'markerfacecolor','r');
  plot(x,pmstd,'g');
  if ~isempty(theorout)
    hold on; plot([x(1),x(end)],[theorout(rv),theorout(rv)],'--k');
    ylim(theorout(rv).*[1-frac(rv),1+frac(rv)]); % axis tight
  else
    axis tight
  end
  
  if ~isempty(ya)
    ylim(ya(rv,:))
  end
  
  xlim([x(1),x(end)])
  xlabel(xlstr)
  ylabel(sprintf('%s',RVstr{rv}));
  
  title({[sprintf('%s',titlestr)],...
    [sprintf('%d runs per sample size',maxMC)]})
  text(x(end),ya(rv,1) + abs(ya(rv,1)-ya(rv,2))/10,...
    sprintf('%s theory = %1.2f',RVstr{rv},theorout(rv)),...
    'HorizontalAlignment','right','fontsize',12);
end
