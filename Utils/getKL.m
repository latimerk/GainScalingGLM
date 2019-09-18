function [KL] = getKL(P,Q)

P = P(:)./sum(P);
Q = Q(:)./sum(Q);

if(sum(Q == 0 & P > 0))
    warning('KL divergence cannot be calculated: Q does not share support with P');
    KL = nan;
    return;
end

T = P > 0;

KL = P(T)'*(log2(P(T)) - log2(Q(T)));