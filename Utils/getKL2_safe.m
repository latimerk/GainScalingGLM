function [KL] = getKL2_safe(P,Q)
% Q(P==0) = 0;
% P(Q==0) = 0;
Q = max(eps,Q);
P = max(eps,P);

KL = 1/2*(getKL(P,Q) + getKL(Q,P));