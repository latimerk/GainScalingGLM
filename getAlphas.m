load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
%load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');
addpath Utils/
spkHistLengths(1) = 10;

StimLevels = [1;1.3;1.6;2];
StimPeriods = [1;2;4;8;16;32;64];

HH = length(spkHistLengths);
saveFigs = false;
NB = 30;


pp = StimPeriods./NB;


alphas.filter.square = nan(HH,3);
alphas.filter.flat   = nan(HH,3);
alphas.filter.sine   = nan(HH,3);
alphas.filter.hh     = nan(1,3);
alphas.filter.hh_def = nan(1,3);

for cc = 1:3
    yy_hh = fits_time.cycleAvg_3AHP(:,:,cc+1);
    gg = fits_time.cycleAvg_Gain(:,:,cc+1);

    alphas.filter.hh(cc) = fitAlphaSteps(yy_hh,gg-1,pp);

    yy_hh = fits_time.cycleAvg_Default(:,:,cc+1);
    alphas.filter.hh_def(cc) = fitAlphaSteps(yy_hh,gg-1,pp);


    for hh = 1:HH
        yy_square = fits_time.cycleAvg_glm_3AHP(:,:,cc+1,hh);
        yy_flat = fits_time.cycleAvg_glm_flatTrain_3AHP(:,:,cc+1,hh);
        yy_sine = fits_time.cycleAvg_glm_sineTrain_3AHP(:,:,cc+1,hh);

        alphas.filter.square(hh,cc) = fitAlphaSteps(yy_square,gg-1,pp);
        alphas.filter.flat(hh,cc)   = fitAlphaSteps(yy_flat,gg-1,pp);
        alphas.filter.sine(hh,cc)   = fitAlphaSteps(yy_sine,gg-1,pp);
    end
end


%%
alphas.gain.square = nan(HH,3);
alphas.gain.flat   = nan(HH,3);
alphas.gain.sine   = nan(HH,3);
alphas.gain.hh     = nan(1,3);
alphas.gain.hh_def = nan(1,3);

alphas.gain.psUsed = 1:7;

X = log(StimPeriods(alphas.gain.psUsed));
X = [-X ones(size(X))];
XX = X'*X;

for cc = 1:3

    bb_hh = XX\(X'*log(fits_phase.gs_3AHP( alphas.gain.psUsed,cc+1)));
    alphas.gain.hh(cc) = bb_hh(1);

    bb_hh = XX\(X'*log(fits_phase.gs_Default( alphas.gain.psUsed,cc+1)));
    alphas.gain.hh_def(cc) = bb_hh(1);

    for hh = 1:HH

        bb_glm = XX\(X'*log(fits_phase.gs_glm_sqTrain_3AHP( alphas.gain.psUsed,cc+1,hh)));
        alphas.gain.square(hh,cc) = bb_glm(1);

        bb_glm = XX\(X'*log(fits_phase.gs_glm_flatTrain_3AHP( alphas.gain.psUsed,cc+1,hh)));
        alphas.gain.flat(hh,cc) = bb_glm(1);

        bb_glm = XX\(X'*log(fits_phase.gs_glm_3AHP( alphas.gain.psUsed,cc+1,hh)));
        alphas.gain.sine(hh,cc) = bb_glm(1);
    end
end


%%
alphas.phase.square = nan(HH,3);
alphas.phase.flat   = nan(HH,3);
alphas.phase.sine   = nan(HH,3);
alphas.phase.hh     = nan(1,3);
alphas.phase.hh_def = nan(1,3);

alphas.phase.psUsed = 1:7;


for cc = 1:3

    gg = fits_phase.phase_Gain(:,cc+1);
    alphas.phase.hh(cc) = mean(2*(fits_phase.phase_3AHP(alphas.phase.psUsed,cc+1) - gg)./pi);
    alphas.phase.hh_def(cc) = mean(2*(fits_phase.phase_Default(alphas.phase.psUsed,cc+1) - gg)./pi);


    for hh = 1:HH

        alphas.phase.square(hh,cc) = mean(2*(fits_phase.phase_glm_sqTrain_3AHP(alphas.phase.psUsed,cc+1,hh)-gg)./pi);

        alphas.phase.flat(hh,cc) = mean(2*(fits_phase.phase_glm_flatTrain_3AHP(alphas.phase.psUsed,cc+1,hh)-gg)./pi);

        alphas.phase.sine(hh,cc) = mean(2*(fits_phase.phase_glm_3AHP(alphas.phase.psUsed,cc+1,hh)-gg)./pi);
    end
end

%%

save('Results/Lundstrom/alphas.mat','alphas','spkHistLengths','StimPeriods','StimLevels');