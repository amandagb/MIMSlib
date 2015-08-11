function [flag,rise_temp,fall_temp] = find_rise(y,line_num,top_thresh,bot_thresh,varargin)
%% Script: align_rise = find_rise(y,line_num,template,avgal)
% Description: Takes the OES standard text file with structure
% LABELline#,MM/DD/YYYY,
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% y:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  29 Mar  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     Seperating subfunction from reshape_OEScorr
%  21 June 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     No changes, deleted extraneous code

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end
y = medfilt1(y,9);
dfig = figure;%get(gcf);%
flag = zeros(1,2);
fall_temp = 0;
rise_temp = 0;
for j = 1:2;
  Ny = length(y);
  tfmid = (y < top_thresh & y > bot_thresh);
  tftop = (y >= top_thresh).*3;
  tfy = tfmid + tftop - ones(1,Ny);
  dtfy = [0,diff(tfy)]; 
%   Suppose you divide the space into three regions: I = points < bottom
%   threshold; II - point between bottom and top thresholds; III - points
%   above top threshold. dtfy can take on the following values: dtfy = 1 ->
%   region I TO region II; dtfy = 2 -> region II TO region III; dtfy = 3 ->
%   region I TO region III; dtfy = -1 -> region II TO region I; dtfy = -2
%   -> region III TO region II; dtfy = -3 -> region III TO region I
  non0dtfy = find(dtfy ~= 0);
  
  plot(y,'marker','o','markersize',2); hold on;
  plot(xlim,top_thresh.*ones(1,2),':k');
  plot(xlim,bot_thresh.*ones(1,2),':k');
  plot(((tfy+1)/3).*max(y),'r.'); hold off;
  
  i = 1;
  if ~isempty(non0dtfy) % NEW LINE, added 10 May 2012
    while i <= length(non0dtfy)-1 && ~flag(j)
      if dtfy(non0dtfy(i)) == 1 && dtfy(non0dtfy(i+1)) == 2
        % Case where the data goes from region I to region II and ends in
        % region III
        flag(j) = non0dtfy(i) - 1;
        if j == 1
          rise_temp = y(flag(j):non0dtfy(i+1)); % points in the group which go from region I to region III
        elseif j == 2
          fall_temp = fliplr(y(flag(j):non0dtfy(i+1))); % points in the group which go from region III to region I
        end
      elseif dtfy(non0dtfy(i)) == 3 && abs(sum(dtfy(non0dtfy(i) + (1:3)))) <= 2
        % This logic statement indicates that points when from region I to
        % region III without populating region II in between
        flag(j) = non0dtfy(i) - 1;
        if j == 1
          rise_temp = y(flag(j)+ (0:3));
        elseif j == 2
          fall_temp = fliplr(y(flag(j)+ (0:3)));
        end
      else
        i = i + 1;
      end
    end
    y = fliplr(y);
  else
    flag(j) = 0;
  end
end
if sum(flag) > 0
  flag(j) = length(y)-flag(j)+1-length(fall_temp)+1;
end
plot(y,'marker','o'); hold on;
plot(xlim,top_thresh.*ones(1,2),':k');
plot(xlim,bot_thresh.*ones(1,2),':k');
plot((1:length(rise_temp))+flag(1),rise_temp,'*g');
plot((1:length(fall_temp))+flag(2),fall_temp,'*g');
title(sprintf('Line %d',line_num));
hold off;
if line_num == 1 %mod(line_num,10) == 1
  pause
end

close(dfig)
end
