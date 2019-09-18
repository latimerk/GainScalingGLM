function [JS] = getJS(P,Q)
P = P(:)./sum(P);
Q = Q(:)./sum(Q);

M = 1/2*(P+Q);

JS = 1/2*(getKL(P,M) + getKL(Q,M));