function [KL] = getAM_safe(P,Q,D)

Q = max(eps,Q);
P = max(eps,P);
D = max(eps,D);
P = P(:)./sum(P);
Q = Q(:)./sum(Q);
D = D(:)./sum(D);

if(sum(Q == 0 & P > 0))
    warning('KL divergence cannot be calculated: Q does not share support with P');
    KL = nan;
    return;
end

T = P > 0;

KL = D(T)'*(log2(P(T)) - log2(Q(T)));