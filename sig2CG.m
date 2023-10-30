function [ sig2 ] = sig2CG( N, alpha)
%
% Type: [ sig2 ] = sig2CG( N, alf)
%
% Inputs:
%
% N   := 1 x 1 window length
% alpha := 1 x 1 parameter of the Confined Gaussian Window
%
% Outputs:
%
% sig2   := 1 x 1 sig^2 of the CG window
%
% This funciton is only used by replacewinbyCG().
%
% Daniel Haegele, Ruhr-Universitaet Bochum, August 2015
  w = Cgausswin(N,alpha);
  norm2 = w'*w;
  pos = (1:N)';
  sig2 = (w'*diag(pos.^2)*w/norm2 - (w'*diag(pos)*w/norm2)^2);
end

