function [KL] = getKL2_nonparam(P,Q,k)

if(nargin < 3)
    KL = 1/2*(getKL_nonparam(P,Q) + getKL_nonparam(Q,P));
else
    KL = 1/2*(getKL_nonparam(P,Q,k) + getKL_nonparam(Q,P,k));
end