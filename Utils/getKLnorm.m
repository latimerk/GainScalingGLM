function [kl] = getKLnorm(m0,m1,sig0,sig1)

addpath ~/gitCode/GLM/code/utils/
kl = trace(sig1\sig0) + (m1-m0)'*(sig1\(m1-m0))-length(m0) + logdet(sig1) - logdet(sig0);

kl = kl./log(2);