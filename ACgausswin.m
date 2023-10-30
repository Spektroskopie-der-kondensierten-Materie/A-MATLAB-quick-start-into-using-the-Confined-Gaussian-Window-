function [w] = ACgausswin(N,alf)
%
% Type: [w] = ACgausswin(N,alf);
%
% Inputs:
%
% N   := 1 x 1 number of points in window
% alf := 1 x 1 bandwidth parameter (default : 2.5)
%              large alpha - temporal width of window is small and
%              spectral bandwidth is large.
%
% Outputs:
%
% w   := n x 1 sampled window
%
%
% Compute Approximate Confined Gaussian (ACG) window - the windows with
% the best possible RMS time-bandwidth product.

% The function ACgausswin() is used by u = replacewinbyAGC(w) to compute
% an AGCwindow with same temporal width and same vector norm, i.e. 
% u'*u will be equal to w'*w.
%
% The implementation of ACgausswin() was derived from gausswin().
% Be aware that the RMS temporal width of the two windows is
% different for the same parameter alpha. gausswin() yields a larger temporal
% width than ACgausswin(), but ACgausswin() still yields a small 
% RMS frequency width, i.e. its time-bandwidth product is better. 
% % For small temporal width, i.e. alf > 7, ACgausswin and gausswin yield
% identical windows.
% The main advantage of ACgausswin() over gausswin() is that ACgausswin()
% yields an almost perfect approximation of the Confined Gaussian window
% which posseses the best possible RMS time-bandwidth product for a given
% RMS temporal width.
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

msg=nargchk(1,2,nargin);
if ~isempty(msg)
 error(msg);
end
if nargin < 2 | isempty(alf),
 alf=2.5;
end
%


G = @(x) exp(-(1/2)*(alf*(x-(N-1)/2)/((N-1)/2)).^2);
p=(0:(N-1))';
if alf < 7
    w= G(p) - G(-1/2)*(G(p+N) + G(p-N))/(G(-1/2+N)+G(1/2-N));
else
    w = G(p);
end    
%%

