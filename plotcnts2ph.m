function Iph = plotcnts2ph(d,imgstr)
%% function: plotcnts2ph
% Description:
% INPUTS ----------------------------------------------------------------
%
% OUTPUTS ---------------------------------------------------------------
% comp
%
%  Date           Author            E-mail                      Version
%  25 Nov  2014   Amanda Balderrama amanda.gaudreau@gmail.com     0
imgstr = {'luminescent'};%'fluorescentreference',...
%   'luminescent',...
%   'photograph',...
%   'readbiasonly'};
clknumind = nonemptystrind(d.SequenceInfo.info(:,2),'ClickNumber ');
nimgs = length(clknumind);
cnt2phC = zeros(1,nimgs);
for i = 1:length(imgstr)
  imtag = imgstr{i};
  dimstruct = getfield(d,imtag);
  dIcnt = getfield(dimstruct,'image');
  dimfields = fieldnames(dimstruct);
  FOVid = nonemptystrind(lower(dimfields),'fieldofview');
  FOVnum = getfield(dimstruct,dimfields{FOVid});
  fnumid = nonemptystrind(lower(dimfields),'fnumber');
  fnumnum = getfield(dimstruct,dimfields{fnumid});
  
  for n = 1:nimgs
    if rem(FOVnum(n),1) ~= 0
      cnt2phCind = nonemptystrind(d.ClickTable(:,2),sprintf('Coef C-ccd at FOV %1.1f, f%d',FOVnum(n),fnumnum(n)));
    else
      cnt2phCind = nonemptystrind(d.ClickTable(:,2),sprintf('Coef C-ccd at FOV %d, f%d',FOVnum(n),fnumnum(n)));
    end
    cnt2phC(n) = str2num(d.ClickTable{cnt2phCind,3});
  end
end
end