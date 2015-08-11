function [dout] = keeplargestCC(dbin)
%% Script: [dout] = keeplargestCC(dbin)
% Description: Keeps the largest connected component from a binary mask
% Example:  
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% dbin: binary image
% OUTPUTS ---------------------------------------------------------------
% dout: largest component of original input image
%
%  Date           Author            E-mail                      Version
%  1  Dec  2014   Amanda Balderrama amanda.gaudreau@gmail.com     1

cc = bwconncomp(dbin);
Lcc = cellfun(@(x) length(x),cc.PixelIdxList);
[leni,keepi] = max(Lcc);
dout = false(size(dbin));
dout(cc.PixelIdxList{keepi}) = true;
end
