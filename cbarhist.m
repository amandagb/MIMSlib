function cbarhist(I,varargin)
%% Script: cbarhist(I,varargin)
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%
% Date           Author            E-mail                      Version
% 04 June 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

figure;
% Determine user options
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('nbins',PropertyNames)
  ncol = PropertyVal{strmatch('nbins',PropertyNames)};
else
  ncol = 64;
end

if strmatch('colormap',PropertyNames)
  colmapstr = PropertyVal{strmatch('colormap',PropertyNames)};
else
  colmapstr = 'gray';
end

if strmatch('CLim',PropertyNames)
  CLim = PropertyVal{strmatch('CLim',PropertyNames)};
else
  CLim = [min(I(:)),max(I(:))];
end

Cmap = custom_colormaps(colmapstr,ncol);

% Prepare image variables
Imin = CLim(1);
Imax = CLim(2);
Irng = CLim(2) - CLim(1);
binwidth = Irng/(ncol);
bincenters = [(CLim(1) + binwidth/2):binwidth:CLim(2)];
range = [CLim(1),CLim(2)];

y = hist(I(:),bincenters);
y = y(:);
bar(bincenters,y,'hist')
hist_axes = gca;
h_fig = ancestor(hist_axes,'figure');
axis tight
limits = axis(hist_axes);
axis(hist_axes,limits);
original_axes_pos = get(hist_axes,'Position');
hist_axes_units_old = get(hist_axes,'units');
set(hist_axes,'Units','Normalized');
% Get axis position and make room for color stripe.
pos = get(hist_axes,'pos');
stripe = 0.075;
set(hist_axes,'pos',[pos(1) pos(2)+stripe*pos(4) pos(3) (1-stripe)*pos(4)])
set(hist_axes,'Units',hist_axes_units_old);
set(hist_axes,'xticklabel','')
stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
  'Position', [pos(1) pos(2) pos(3) stripe*pos(4)]);
limits = axis(stripe_axes);
limits(1:2) = range;
xdata = [bincenters(1),bincenters(end)];%[CLim(1)+binwidth/2,CLim(2)];
C = (1:ncol)/ncol;
Cm = reshape(Cmap,[1,ncol,3]);%Cm = repmat(C, [1 1 3]);
imagesc(xdata,[0 1],Cm,'Parent',stripe_axes);
set(stripe_axes,'yticklabel','')
axis(stripe_axes,limits);

% Put a border around the stripe.
line(limits([1 2 2 1 1]),limits([3 3 4 4 3]),...
  'LineStyle','-',...
  'Parent',stripe_axes,...
  'Color',get(stripe_axes,'XColor'));
set(h_fig,'CurrentAxes',hist_axes);

% Tag for testing.
set(stripe_axes,'tag','colorstripe');

wireHistogramAxesListeners(hist_axes,stripe_axes,original_axes_pos);

% Link the XLim of histogram and color stripe axes together.
% In calls to imhist in a tight loop, the histogram and colorstripe axes
% are destroyed and recreated repetitively. Use linkprop rather than
% linkaxes to link xlimits together to solve deletion timing problems.
h_link = linkprop([hist_axes,stripe_axes],'XLim');
setappdata(stripe_axes,'linkColorStripe',h_link);

%%%
%%% Function wireHistogramAxesListeners
%%%
  function wireHistogramAxesListeners(hist_axes,stripe_axes,original_axes_pos)
    
    % If the histogram axes is deleted, delete the color stripe associated with
    % the histogram axes.
    cb_fun = @(obj,evt) removeColorStripeAxes(stripe_axes);
    lis.histogramAxesDeletedListener = iptui.iptaddlistener(hist_axes,...
      'ObjectBeingDestroyed',cb_fun);
    
    % This is a dummy hg object used to listen for when the histogram axes is cleared.
    deleteProxy = text('Parent',hist_axes,...
      'Visible','Off', ...
      'Tag','axes cleared proxy',...
      'HandleVisibility','off');
    
    % deleteProxy is an invisible text object that is parented to the histogram
    % axes.  If the ObjectBeingDestroyed listener fires, the histogram axes has
    % been cleared. This listener is triggered by newplot when newplot clears
    % the current axes to make way for new hg objects being drawn. This
    % listener does NOT fire as a result of the parent axes being deleted.
    prox_del_cb = @(obj,evt) histogramAxesCleared(obj,stripe_axes,original_axes_pos);
    lis.proxydeleted = iptui.iptaddlistener(deleteProxy,...
      'ObjectBeingDestroyed',prox_del_cb);
    
    setappdata(stripe_axes,'ColorStripeListeners',lis);
  end

%%%
%%% Function removeColorStripeAxes
%%%
  function removeColorStripeAxes(stripe_axes)
    
    if ishghandle(stripe_axes)
      delete(stripe_axes);
    end
  end

%%%
%%% Function histogramAxesCleared
%%%
  function histogramAxesCleared(hDeleteProxy,stripe_axes,original_axes_pos)
    
    removeColorStripeAxes(stripe_axes);
    
    h_hist_ax = get(hDeleteProxy,'parent');
    set(h_hist_ax,'Position',original_axes_pos);
  end
end