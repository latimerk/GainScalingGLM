function [kl] = getKL_nonparam(P,Q,k)
%Implementation of k-nn estimator of KL divergence with dimension=1 from:
%
%Estimation of Information Theoretic Measures for Continuous Random Variables
%Fernando Perez-Cruz, NIPS 2009
%https://papers.nips.cc/paper/3417-estimation-of-information-theoretic-measures-for-continuous-random-variables.pdf
%
% code by Kenneth Latimer
n = length(P);
m = length(Q);

if(nargin < 3)
    k = min([floor(n/2) floor(m/2) 20]);
end
   
if(n - 1 < k || m < k) 
    error('Too few observations for given k');
end

if(n < 10 || m < 10)
    warning('Low number of observations!');
end

[~,r_k] = knnsearch(P,P,'K',k+1);
% r_k = sum(r_k(:,2:end),2);
r_k = r_k(:,end);

[~,s_k] = knnsearch(Q,P,'K',k);
% s_k = sum(s_k(:,1:end),2);
s_k = s_k(:,end);

kl = sum(log(s_k) - log(r_k))/n + log(m/(n-1));