function [ u, alpha] = replacewinbyCG( w )
%
% Type: [u,alf] = replacewinbyACG( w);
%
% Inputs:
%
% w   := n x 1 sampled window
%
% Outputs:
%
% u   := n x 1 sampled window
%
% The function u = replacewinbyAGC(w) uses ACgausswin() to compute
% an AGCwindow with same temporal width and same vector norm, i.e. 
% u'*u is equal to w'*w.
%
% The Confined Gaussian window and the ACG window was introduced 
% in the paper
%   Sebatian Starosielec and Daniel Haegele
%     Discrete-time windows with minimal RMS bandwidth
%     for given RMS temporal width, 
%   Signal Processing 102, 240 (2014) 
% 
% Please cite this paper in your work when you use the Confined or
% Approximate Confined Gaussian Window
%
% Daniel Haegele, Ruhr-Universitaet Bochum, August 2015
  N = length(w);
  norm2 = w'*w;
  pos = (1:N)';
  % calculate the temporal sig^2 of the input window 
  sig2 = (w'*diag(pos.^2)*w/norm2 - (w'*diag(pos)*w/norm2)^2);
  if sqrt(sig2)/N < 0.5
  % alf0 = (N-1)/(2*sqrt(2*sig2));
  alpha0 = 0.4;
   % sig2ACG(N,alf) computes the temporal sig^2 of the ACG window for
    % window length N and window parameter alpha (variable alf).
    x = fminsearch(@(xx) ((sig2CG(N,xx))-sig2)^2,alpha0,optimset('TolX',1e-10));
    uh = Cgausswin(N,x);
    u = uh/sqrt(uh'*uh)*sqrt(norm2);
    alpha = x;
  else  
   u = w ;
   alpha = -1.0;   
  end
end  

