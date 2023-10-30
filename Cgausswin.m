function [ w ] = Cgausswin( N, alphaInUnits )
%
% Type: [w] = Cgausswin(N,alf);
%
% Inputs:
%
% N   := 1 x 1 number of points in window
% alphaInUnits := 1 x 1 Is the paramater alpha in units of alphabar 
% as defined in the original definition by Starosielec and Haegele in
% reference [1]
% Outputs:
%
% w   := n x 1 sampled window
%
% The Confined Gaussian window and the Approximate Confined Gaussian window
% was introduced in the paper

% [1]  Sebatian Starosielec and Daniel Haegele
%     Discrete-time windows with minimal RMS bandwidth
%     for given RMS temporal width, 
%   Signal Processing 102, 240 (2014) 
% 
% Please cite this paper in your work when you use the Confined or
% Approximate Confined Gaussian Window
    abar = (10/N)^4/4; %alphabar
    alpha = alphaInUnits*abar;
    % Calulate the Confined Gaussian Window along the procedure given in
    % [1]. 
    % Calculation of T and P as defined in [1]
    T = zeros(N,N);
    P = zeros(N,N);
    for k=1:N
        T(k,k) = (k - (N+1)/2)^2;
        for l=1:N
            if k ~= l
                P(k,l) = 2*(-1)^(k-l)/(k-l)^2;
            else
                P(k,l) = pi^2/3;
            end
        end
    end
    opts.maxit = 10000;
    % Calculate the eigenvector with the lowest eigenvalue
    [w,lambda] = eigs(P + alpha*T, 1, 'sa', opts); 
    w = abs(w);
end

