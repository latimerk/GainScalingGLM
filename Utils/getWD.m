function [EMD] = getWD(P,Q,xx)

addpath ~/gitCode/GainScaling/Utils/EMD/

P = P(:)'./sum(P);
Q = Q(:)'./sum(Q);

xx = xx(:);

C = abs(xx-xx');

EMD = emd_mex(P,Q,C); 