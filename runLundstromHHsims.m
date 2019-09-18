%Simulate HH neurons as in
% Lundstrom, B. N., Higgs, M. H., Spain, W. J., and Fairhall, A. L. (2008).  Fractional differentiation by neocortical pyramidal neurons.
%   Nature neuroscience
%

rng(08022018);



debug = false;




sigMuRatio = 4.0;
StimLevels = [1.0; 1.3; 1.6; 2.0];
StimPeriods = [1; 2; 4; 8; 16; 32; 64];
targetRate = 10.0;
downsampleRate = 100;


if(debug)
    BUFFER = 1;
    L = 100000-2;
    L_init = 100000-2;
else
    BUFFER = 60*1000*downsampleRate;
    L = 1000*(max(StimPeriods)*50)*downsampleRate;
    L_init = 100*1000*downsampleRate;
end
h = 1.0/(downsampleRate);
T_range = BUFFER +[1 L];
T_range_init = BUFFER +[1 L_init];
xs = 0.0:0.01:5.5;
frFunctions = nan(length(xs),3);


T_range_downsampled = (ceil(T_range./downsampleRate));

X1 = randn((ceil((L+2*BUFFER)/downsampleRate)),1);


y = zeros(length(X1)*downsampleRate,1);
for ii = 1:length(X1)
    y((1:downsampleRate) + (ii-1)*downsampleRate) = X1(ii);
end

y_init = y(1:(T_range_init(end) + 1));


NP = length(StimPeriods);
NL = length(StimLevels);




%%
simTypes = {"Lundstrom3AHP","Lundstrom2Na","LundstromDefault"};
fprintf("Getting baseline rates...\n");
parfor ii = 1:length(simTypes)
    [~,frFunctions(:,ii),rs(ii)] = getStimulusRate(simTypes{ii},xs,targetRate,sigMuRatio,T_range_init,y_init,h);
end
stim_mu_3AHP = rs(1);
stim_mu_2Na = rs(2);
stim_mu_Default = rs(3);



fprintf("Starting HH sims.\n")

sqGainFunc   = @(t,Period,Height) round(sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
sineGainFunc = @(t,Period,Height)      (sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;



sts_0_sq_2Na = cell(NP,NL);
sts_0_sine_2Na = cell(NP,NL);
sts_0_sq_3AHP = cell(NP,NL);
sts_0_sine_3AHP = cell(NP,NL);
sts_0_sq_Default = cell(NP,NL);
sts_0_sine_Default = cell(NP,NL);

parfor ii = 1:NP
    sts_0_sq_2Na(ii,:)   = simLundstrom_2Na(y,h,stim_mu_2Na*sigMuRatio,stim_mu_2Na,@(t)sqGainFunc(t,StimPeriods(ii),StimLevels));
    sts_0_sine_2Na(ii,:) = simLundstrom_2Na(y,h,stim_mu_2Na*sigMuRatio,stim_mu_2Na,@(t)sineGainFunc(t,StimPeriods(ii),StimLevels));
    
    sts_0_sq_3AHP(ii,:)   = simLundstrom_3AHP(y,h,stim_mu_3AHP*sigMuRatio,stim_mu_3AHP,@(t)sqGainFunc(t,StimPeriods(ii),StimLevels));
    sts_0_sine_3AHP(ii,:) = simLundstrom_3AHP(y,h,stim_mu_3AHP*sigMuRatio,stim_mu_3AHP,@(t)sineGainFunc(t,StimPeriods(ii),StimLevels));
    
    sts_0_sq_Default(ii,:)   = simLundstrom_Default(y,h,stim_mu_Default*sigMuRatio,stim_mu_Default,@(t)sqGainFunc(t,StimPeriods(ii),StimLevels));
    sts_0_sine_Default(ii,:) = simLundstrom_Default(y,h,stim_mu_Default*sigMuRatio,stim_mu_Default,@(t)sineGainFunc(t,StimPeriods(ii),StimLevels));
end


sts_sq_2Na = cell(NP,NL);
sts_sine_2Na = cell(NP,NL);
sts_sq_3AHP = cell(NP,NL);
sts_sine_3AHP = cell(NP,NL);
sts_sq_Default = cell(NP,NL);
sts_sine_Default = cell(NP,NL);

for ii = 1:NP
    for jj = 1:NL
        sts_sq_2Na{ii,jj}       = ceil(sts_0_sq_2Na{ii,jj}      ./downsampleRate);
        sts_sine_2Na{ii,jj}     = ceil(sts_0_sine_2Na{ii,jj}    ./downsampleRate);
        sts_sq_3AHP{ii,jj}      = ceil(sts_0_sq_3AHP{ii,jj}     ./downsampleRate);
        sts_sine_3AHP{ii,jj}    = ceil(sts_0_sine_3AHP{ii,jj}   ./downsampleRate);
        sts_sq_Default{ii,jj}   = ceil(sts_0_sq_Default{ii,jj}  ./downsampleRate);
        sts_sine_Default{ii,jj} = ceil(sts_0_sine_Default{ii,jj}./downsampleRate);
    end
end

save("Results/Lundstrom/LundstromSims.mat","-v7.3","X1","StimLevels","StimPeriods","T_range_downsampled","sigMuRatio","stim_mu_2Na",...
                                 "stim_mu_3AHP","stim_mu_Default","targetRate","frFunctions",...
                                 "sts_sq_2Na","sts_sq_3AHP","sts_sq_Default",...
                                 "sts_sine_2Na","sts_sine_3AHP","sts_sine_Default","sqGainFunc","sineGainFunc","debug","h","downsampleRate","BUFFER");

