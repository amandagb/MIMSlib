function new_hist = build_rgbpdf(rgb_colvec,nbins)
%% build_rgbpdf
% Description:
% INPUTS ----------------------------------------------------------------
%
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author              E-mail                      Version
%  26 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%     

if ~exist('nbins');
  nbins = 256;
end

scalefact = 256/nbins;
rgb_colvec = rgb_colvec./scalefact;
rgb_colvec = floor(rgb_colvec);
[Urgb,ia,ic] = unique(rgb_colvec,'rows');
nrgb = hist(ic,1:max(ic));
indhist = sub2ind([nbins,nbins,nbins],Urgb(:,1)+1,Urgb(:,2)+1,Urgb(:,3)+1);
new_hist = uint32(zeros([nbins,nbins,nbins]));
new_hist(indhist) = uint32(nrgb(:));

end