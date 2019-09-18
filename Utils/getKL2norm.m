function [kl] = getKL2norm(m0,m1,sig0,sig1)

kl = 1/2*(getKLnorm(m0,m1,sig0,sig1) + getKLnorm(m1,m0,sig1,sig0));