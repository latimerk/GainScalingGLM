function [KL] = getKL2(P,Q)

KL = 1/2*(getKL(P,Q) + getKL(Q,P));