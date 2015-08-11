function acceptpoint = track_acceptancesa(optimValues,newx,newfval)
%ACCEPTANCESA Acceptance function for simulated annealing solver
%   ACCEPTPOINT = ACCEPTANCESA(optimValues,newX,newfval) uses the
%   change in function values between the current point and new point to
%   determine whether the new point is accepted or not.
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration
%             t0: start time
%              k: annealing parameter
%
%   NEWX: new point
%
%   NEWFVAL: function value at NEWX
%
%   Example:
%    Create an options structure using ACCEPTANCESA as the annealing
%    function
%    options = saoptimset('AcceptanceFcn',@acceptancesa);

%   Copyright 2006-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/05/10 17:16:41 $
%   17 Jan  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%   12 May  2015   Amanda Gaudreau   amanda.gaudreau@gmail.com     2

cdir = loaddirfun;
save2fldr = cdir.sapath;
acceptall = 0;
useh = 'h1';
dateinfo = cdir.dateyyyymmddstr;
seededuni = 0;
expnum = 1;
nvar = numel(optimValues.x);

delE = newfval - optimValues.fval;

% If the new point is better accept it
if delE < 0
  acceptpoint = true;
  h1 = NaN; h2 = NaN; num_uni = NaN;
  % Otherwise, accept it randomly based on a Boltzmann probability density
else
  h1 = 1/(1+exp(delE/max(optimValues.temperature)));
  h2 = exp(-delE/max(optimValues.temperature));
  switch useh
    case 'h1'
      h = h1;%1/(1+exp(delE/max(optimValues.temperature)));
    case 'h2'
      h = h2;%exp(-delE/max(optimValues.temperature));
  end
  if seededuni
    %f =  'C:\Users\Amanda\Documents\Dropbox\MADLab Research\Data\simanneal Experiments\seed0_5kby1_rand.csv';
    f = strcat(cdir.sapath,'seed0_5kby1_rand.csv');
    num_uni = csvread(f,optimValues.iteration,0,[optimValues.iteration,0,optimValues.iteration,0])';
  else
    num_uni = rand;
  end
  if acceptall
    acceptpoint = true;
  elseif h > num_uni
    acceptpoint = true;
  else
    acceptpoint = false;
  end
end

acceptfilename = sprintf('%s_nvar%dPaccept%s_exp%d.csv',dateinfo,nvar,useh,expnum);
fname = sprintf('%s%s',save2fldr,acceptfilename);
dlmwrite(fname,[optimValues.fval, optimValues.x(:)', newfval, newx(:)', delE, ...
  h1, num_uni, max(optimValues.temperature), acceptpoint, datenum(now)],'-append');

% -------------------------------------------------------------------------
% COMMENTED OUT 9 MAY 2013 --- this file is not used and can be derrived
%   from the Paccept file
% -------------------------------------------------------------------------
% bestfilename = sprintf('%s_nvar%dbest%s_exp%d.csv',dateinfo,nvar,useh,expnum);
% % bestfilename = sprintf('%s_nvar%d_best%s_exp%d.csv',dateinfo,nvar,useh,expnum);
% bestfname = sprintf('%s%s',fldr,bestfilename);
% if newfval <= optimValues.bestfval
%     dlmwrite(bestfname,[optimValues.iteration, optimValues.temperature(1),...
%         newfval,newx],'-append');
% end
