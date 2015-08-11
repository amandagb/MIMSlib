function strind = nonemptystrind(str,pattern)
%% function: nonemptystrind
% Description: 
% INPUTS ----------------------------------------------------------------
%
% OUTPUTS ---------------------------------------------------------------
% comp    
%
%  Date           Author            E-mail                      Version
%  24 Nov  2014   Amanda Balderrama amanda.gaudreau@gmail.com     0

strind = find(cellfun(@(x) ~isempty(x),strfind(str,pattern)));

end