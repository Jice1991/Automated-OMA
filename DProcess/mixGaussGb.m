function [label, Theta, w, llh, cluster] = mixGaussGb(X, opt)
% Collapsed Gibbs sampling for Dirichlet process (infinite) Gaussian mixture model (a.k.a. DPGM). 
% This is a wrapper function which calls underlying Dirichlet process mixture model.
% Input: 
%   X: d x n data matrix
%   opt(optional): prior parameters
% Output:
%   label: 1 x n cluster label
%   Theta: 1 x k structure of trained Gaussian components
%   w: 1 x k component weight vector
%   llh: loglikelihood

[d,n] = size(X);
mu = mean(X,2);
Xo = bsxfun(@minus,X,mu);
s = sum(Xo(:).^2)/(d*n);
if nargin == 1
    kappa0 = 0.01; % or smaller
    m0 = mean(X,2);
    nu0 = 50; % original is d
    S0 = s*eye(d);
    alpha0 = 50; % concentrating factor in GP very important!
else
    kappa0 = opt.kappa;
    m0 = opt.m;
    nu0 = opt.nu;
    S0 = opt.S;
    alpha0 = opt.alpha;
end
prior = GaussWishart(kappa0,m0,nu0,S0); % base measure, conjugate prior
[label, Theta, w, llh, cluster] = mixDpGb(X,alpha0,prior);


