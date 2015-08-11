function linespersamp = find_linespersamp(data);
%% Script: linespersamp = find_linespersamp(data);
% Description: Determines the number of lines per sample from raw OES data
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% data: M x N matrix of the template element (usualy C_)
% varargin - 'PropertyName','PropertyValue'
% 
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  21 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  20 June 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2

meanI = mean(data(:));
stdI = std(data(:));
data = mean(data,1);
data = medfilt1(data,5);
srtdata = sort(data);
srtdata = unique(srtdata);
[v,srti] = min(abs((meanI-stdI/2)-srtdata));
T = srtdata(srti);
lenCC = 0;
mindiffCC = -1;
while sum(abs(lenCC-mean(lenCC)) > mindiffCC)
  indgtT = data > T;
  CC = bwconncomp(indgtT);
  lenCC = cellfun(@(x) length(x),CC.PixelIdxList);
  mindiffCC = min(lenCC)*0.2;%pntsperobj - pntsperobj*0.3;
  if sum(abs(lenCC-mean(lenCC)) > mindiffCC)
    srti = srti - 1;
    T = srtdata(srti);
    disp('Threshold needs to be adjusted')
  end
end
linespersamp = CC.NumObjects;
end

