load('Results/Mease/MeaseSims.mat','X1','downsampleRate','sigMuRatio','StimLevels','sigMuRatio','stim_mu');


h = 1.0/(downsampleRate);
Gs = 600:100:2000;
nn = 5;
kk = 5;

G_Na = Gs(nn);
G_K  = Gs(kk);



X1 = X1(1:(12*1000));

y = zeros(length(X1)*downsampleRate,1);
for ii = 1:length(X1)
    y((1:downsampleRate) + (ii-1)*downsampleRate) = X1(ii);
end

stim_dc = stim_mu(nn,kk);

sigs = stim_mu(nn,kk)*StimLevels*sigMuRatio; 
[sts,V] = simMease(y,h,sigs(:),stim_dc,G_Na,G_K); 
save('Results/Mease/HHexample.mat','G_K','G_Na','downsampleRate','StimLevels','y','V','sts','stim_dc','sigs');