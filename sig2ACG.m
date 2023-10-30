function [ sig2 ] = sig2ACG( N, alf)
%
% Type: [ sig2 ] = sig2ACG( N, alf)
%
%
% Inputs:
%
% N   := 1 x 1 window length
% alf := 1 x 1 parameter of the Approximate Confined Gaussian Window
%
% Outputs:
%
% sig2   := 1 x 1 sig^2 of the ACG window
%
% This funciton is only used by replacewinbyACG
%
% Daniel Haegele, Ruhr-Universitaet Bochum, August 2015
  w = ACgausswin(N,alf);
  norm2 = w'*w;
  pos = (1:N)';
  sig2 = (w'*diag(pos.^2)*w/norm2 - (w'*diag(pos)*w/norm2)^2);
end

