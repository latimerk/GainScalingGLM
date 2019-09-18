load('Results/Lundstrom/LundstromSims.mat','X1','downsampleRate','sigMuRatio','StimLevels','sigMuRatio','stim_mu_3AHP','stim_mu_Default');


h = 1.0/(downsampleRate);

StimPeriod = 1;
StimLevel = StimLevels(end);



BUFFER = 2*1000*downsampleRate;



sqGainFunc   = @(t,Period,Height) round(sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
sineGainFunc = @(t,Period,Height)      (sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;




X1 = X1(1:(10*1000));

y = zeros(length(X1)*downsampleRate,1);
for ii = 1:length(X1)
    y((1:downsampleRate) + (ii-1)*downsampleRate) = X1(ii);
end




[sts_sq  ,V_sq  ] = simLundstrom_3AHP(y,h,stim_mu_3AHP*sigMuRatio,stim_mu_3AHP,@(t)sqGainFunc(t,StimPeriod,StimLevel));
[sts_sine,V_sine] = simLundstrom_3AHP(y,h,stim_mu_3AHP*sigMuRatio,stim_mu_3AHP,@(t)sineGainFunc(t,StimPeriod,StimLevel));

y_sine = y.*sineGainFunc((1:length(y))',StimPeriod,StimLevel)*stim_mu_3AHP*sigMuRatio + stim_mu_3AHP;
y_sq   = y.*sqGainFunc(  (1:length(y))',StimPeriod,StimLevel)*stim_mu_3AHP*sigMuRatio + stim_mu_3AHP;

    
save('Results/Lundstrom/HHexample.mat','downsampleRate','StimLevel','V_sq','V_sine','sts_sq','sts_sine','StimPeriod','y_sq','y_sine');