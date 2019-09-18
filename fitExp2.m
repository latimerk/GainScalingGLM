function [tau,gain,dc,ff_fit] = fitExp2(ff,dt)

ff(isnan(ff)) = 0;

tts = (0:(length(ff)-1))*dt;

tau_0 = 1./prctile(tts,25);
gain_0  = ((ff(1)-ff(end)));
dc    = mean(ff(end-2:end));

w_0 = [tau_0;gain_0];


opts = optimoptions('fminunc','gradobj','off','display','off','maxiter',2000);

errF = @(ww)sqErrTau(ww,ff-dc,tts);
if(isnan(errF(w_0)))
    fprintf('here\n');
end

w_fit = fminunc(errF,w_0,opts);

tau = 1./w_fit(1);
gain = w_fit(2);
ff_fit = gain*exp(-tts./tau) + dc;

end

function [f] = sqErrTau(w,ff,tts)


f_fit = w(2)*exp(-tts.*w(1));
f = sum((ff(:)-f_fit(:)).^2);
end
