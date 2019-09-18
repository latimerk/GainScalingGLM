function [alpha,gain,dc,ff_fit] = fitAlphaSteps(ff,gg,hh)



alpha_0 = 0.25;

w_0 = log(alpha_0);


opts = optimoptions('fminunc','gradobj','off','display','off','maxiter',2000);

errF = @(ww)sqErrFD(ww,ff,gg,hh);
if(isnan(errF(w_0)))
    fprintf('here\n');
end

w_fit = fminunc(errF,w_0,opts);

alpha = exp(w_fit);
[~,ff_fit,gain,dc] = errF(w_fit);
end

function [f,f_fit,gain,dc] = sqErrFD(w,ff,gg,hh)

N = length(hh);
alpha = exp(w);
fd = gg;
for ii = 1:N
    fd(:,ii) = fgl_deriv(alpha,gg(:,ii),hh(ii));
end

xx = [fd(:) ones(numel(fd),1)];
bs = (xx'*xx)\(xx'*ff(:));
gain = bs(1);
dc = bs(2);

f_fit = gain*fd + dc;
f = sum((ff(:)-f_fit(:)).^2);
end
