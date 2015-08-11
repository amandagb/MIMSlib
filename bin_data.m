function [Mind,binranges] = bin_data(M,nbins,varargin)
%% Script: [Mind,binranges] = bin_data(M,nbins,varargin)
% Description: Transforms the values in M into bin numbers. The default
% mapping is linear, though user can specify other non-linear mappings (see
% descriptions for inputs).
% Example: 
% Required Functions: 
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% M:  a matrix or vector of real values
% nbins: number of bins (indices)
% varargin - 'PropertyName','PropertyValue'
%   'norm': binary indicating whether to normalize (1) or not (0) the data
%   to the range 0 to 1 [DEFAULT = 1]
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% Mind: indexed matrix or vector where each number corresponds to the
%   associated bin number
% binranges: Gives the maximum value of the original M values which
%   correspond to that bin. For Mind = i corresponding indices in M have
%   the property that M <= binranges(i)
% 
% Date           Author            E-mail                      Version
% 09 May  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

if ~exist('thresh'); thresh = []; end;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('norm',PropertyNames)
  n = PropertyVal{strmatch('norm',PropertyNames)};
else 
  n = 1;
end

Mrng = max(M(:)) - min(M(:));
bin_rngM = Mrng/(nbins-1); % equally spaced bins which cover range of M
Mind = ceil((M-min(M(:)))./bin_rngM + 1/2); % index image which reference which bin the pixel belongs to
binmin_Mval = ([1:50]-3/2)*bin_rngM + min(M(:)); 
binmid_Mval = ([1:50]-2/2)*bin_rngM + min(M(:));
binmax_Mval = ([1:50]-1/2)*bin_rngM + min(M(:));
binranges = binmax_Mval;

end