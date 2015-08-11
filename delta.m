% Chan-Vese Implementation
% By Rami C., Technion - Israel Institute of Technology
% email: rc@tx.technion.ac.il, website: tx.technion.ac.il/~rc
% For academic use only. Please don't distribute this code without a
% permission from the author.
% read: http://tx.technion.ac.il/~rc/Copyright.htm

% regularized delta function
% inputs:
%    phi - current phi
%    epsilon
%
% output:
%    res - result of delta function


function res=delta(phi,epsilon)


res=(1/pi)*epsilon./(epsilon^2+phi.^2);