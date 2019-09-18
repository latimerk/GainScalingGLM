if(~exist('Results/Lundstrom/sts_glm_sq_2Na','var'))
    %load('Results/Lundstrom/LundstromGLMs.mat','sts_glm_sq_2Na','sts_glm_sq_3AHP','sts_glm_sq_Default','glm_ll','glm_h_spk','glm_k_stim','glm_b','simTypes','spkHistBases','spkHistBases_0','spkHistLengths','stimBasis');
%     load('Results/Lundstrom/LundstromGLMs_part1.mat')
    load('Results/Lundstrom/LundstromGLMs_part2.mat')
%     load('Results/Lundstrom/LundstromGLMs_part3.mat')
%     load('Results/Lundstrom/LundstromGLMs_part10.mat')
    load('Results/Lundstrom/LundstromGLMs_part11.mat')
%     load('Results/Lundstrom/LundstromGLMs_part12.mat')
    load('Results/Lundstrom/LundstromGLMs_meta.mat')
end
if(~exist('sts_sine_2Na','var'))
    load('Results/Lundstrom/LundstromSims.mat');
end
%%

NB = 30;

NP = length(StimPeriods);
NL = length(StimLevels);
NH = length(spkHistLengths);

TT = T_range_downsampled(2)-T_range_downsampled(1)+1;

fits_time = struct();

% fits_time.cycleAvg_2Na     = nan(NB,NP,NL);
fits_time.cycleAvg_3AHP    = nan(NB,NP,NL);
fits_time.cycleAvg_Default = nan(NB,NP,NL);
% fits_time.cycleVar_2Na     = nan(NB,NP,NL);
fits_time.cycleVar_3AHP    = nan(NB,NP,NL);
fits_time.cycleVar_Default = nan(NB,NP,NL);
fits_time.cycleAvg_Gain    = nan(NB,NP,NL);

% fits_time.cycleAvg_glm_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_glm_3AHP    = nan(NB,NP,NL,NH);
fits_time.cycleAvg_glm_Default = nan(NB,NP,NL,NH);


% fits_time.cycleAvg_glm_sineTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_glm_sineTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_glm_sineTrain_Default = nan(NB,NP,NL,NH);

% fits_time.cycleAvg_glm_flatTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_glm_flatTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_glm_flatTrain_Default = nan(NB,NP,NL,NH);

% fits_time.cycleVar_glm_sineTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleVar_glm_sineTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleVar_glm_sineTrain_Default = nan(NB,NP,NL,NH);

% fits_time.cycleVar_glm_flatTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleVar_glm_flatTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleVar_glm_flatTrain_Default = nan(NB,NP,NL,NH);


% fits_time.cycleAvg_hspk_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_hspk_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_hspk_Default = nan(NB,NP,NL,NH);

% fits_time.cycleAvg_kstim_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_kstim_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_kstim_Default = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_kstimv_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_kstimv_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_kstimv_Default = nan(NB,NP,NL,NH);




% fits_time.cycleAvg_hspk2_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleAvg_hspk2_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleAvg_hspk2_Default = nan(NB,NP,NL,NH);

% fits_time.cycleVar_hspk_2Na     = nan(NB,NP,NL,NH);
fits_time.cycleVar_hspk_3AHP    = nan(NB,NP,NL,NH);
% fits_time.cycleVar_hspk_Default = nan(NB,NP,NL,NH);



% fits_time.fitsAvg_2Na     = nan(NB,NP,NL);
fits_time.fitsAvg_3AHP    = nan(NB,NP,NL);
fits_time.fitsAvg_Default = nan(NB,NP,NL);
fits_time.fitsAvg_Gain    = nan(NB,NP,NL);

% fits_time.fitsAvg_glm_2Na     = nan(NB,NP,NL,NH);
fits_time.fitsAvg_glm_3AHP    = nan(NB,NP,NL,NH);
% fits_time.fitsAvg_glm_Default = nan(NB,NP,NL,NH);


% fits_time.fitsAvg_glm_sineTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.fitsAvg_glm_sineTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.fitsAvg_glm_sineTrain_Default = nan(NB,NP,NL,NH);
% fits_time.fitsAvg_glm_flatTrain_2Na     = nan(NB,NP,NL,NH);
fits_time.fitsAvg_glm_flatTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_time.fitsAvg_glm_flatTrain_Default = nan(NB,NP,NL,NH);



% fits_time.fitsAvg_hspk_2Na     = nan(NB,NP,NL,NH);
fits_time.fitsAvg_hspk_3AHP    = nan(NB,NP,NL,NH);
% fits_time.fitsAvg_hspk_Default = nan(NB,NP,NL,NH);

% fits_time.tau_2Na     = nan(NP,NL,2);
fits_time.tau_3AHP    = nan(NP,NL,2);
fits_time.tau_Default = nan(NP,NL,2);
fits_time.tau_Gain    = nan(NP,NL,2);


% fits_time.gs_2Na     = nan(NP,NL,2);
fits_time.gs_3AHP    = nan(NP,NL,2);
fits_time.gs_Default = nan(NP,NL,2);
fits_time.gs_Gain    = nan(NP,NL,2);

% fits_time.dc_2Na     = nan(NP,NL,2);
fits_time.dc_3AHP    = nan(NP,NL,2);
fits_time.dc_Default = nan(NP,NL,2);
fits_time.dc_Gain    = nan(NP,NL,2);


% fits_time.tau_glm_2Na     = nan(NP,NL,NH,2);
fits_time.tau_glm_3AHP    = nan(NP,NL,NH,2);
fits_time.tau_glm_Default = nan(NP,NL,NH,2);

% fits_time.gs_glm_2Na     = nan(NP,NL,NH,2);
fits_time.gs_glm_3AHP    = nan(NP,NL,NH,2);
fits_time.gs_glm_Default = nan(NP,NL,NH,2);

% fits_time.dc_glm_2Na     = nan(NP,NL,NH,2);
fits_time.dc_glm_3AHP    = nan(NP,NL,NH,2);
fits_time.dc_glm_Default = nan(NP,NL,NH,2);




% fits_time.tau_glm_sineTrain_2Na     = nan(NP,NL,NH,2);
fits_time.tau_glm_sineTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.tau_glm_sineTrain_Default = nan(NP,NL,NH,2);
% fits_time.tau_glm_flatTrain_2Na     = nan(NP,NL,NH,2);
fits_time.tau_glm_flatTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.tau_glm_flatTrain_Default = nan(NP,NL,NH,2);

% fits_time.gs_glm_sineTrain_2Na     = nan(NP,NL,NH,2);
fits_time.gs_glm_sineTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.gs_glm_sineTrain_Default = nan(NP,NL,NH,2);
% fits_time.gs_glm_flatTrain_2Na     = nan(NP,NL,NH,2);
fits_time.gs_glm_flatTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.gs_glm_flatTrain_Default = nan(NP,NL,NH,2);

% fits_time.dc_glm_sineTrain_2Na     = nan(NP,NL,NH,2);
fits_time.dc_glm_sineTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.dc_glm_sineTrain_Default = nan(NP,NL,NH,2);
% fits_time.dc_glm_flatTrain_2Na     = nan(NP,NL,NH,2);
fits_time.dc_glm_flatTrain_3AHP    = nan(NP,NL,NH,2);
% fits_time.dc_glm_flatTrain_Default = nan(NP,NL,NH,2);



% fits_time.tau_hspk_2Na     = nan(NP,NL,NH,2);
fits_time.tau_hspk_3AHP    = nan(NP,NL,NH,2);
% fits_time.tau_hspk_Default = nan(NP,NL,NH,2);

% fits_time.gs_hspk_2Na     = nan(NP,NL,NH,2);
fits_time.gs_hspk_3AHP    = nan(NP,NL,NH,2);
% fits_time.gs_hspk_Default = nan(NP,NL,NH,2);

% fits_time.dc_hspk_2Na     = nan(NP,NL,NH,2);
fits_time.dc_hspk_3AHP    = nan(NP,NL,NH,2);
% fits_time.dc_hspk_Default = nan(NP,NL,NH,2);

xx_isi = 0:1:1000;
T_isi = length(xx_isi);
fits_time.ISIs.xx = xx_isi;
fits_time.ISIs.hh_3AHP    = nan(T_isi,NP,NL);
% fits_time.ISIs.hh_2Na     = nan(T_isi,NP,NL);
fits_time.ISIs.hh_Default = nan(T_isi,NP,NL);
fits_time.ISIs.glm_3AHP   = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_2Na    = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_Default = nan(T_isi,NP,NL,NH);
fits_time.ISIs.glm_sineTrain_3AHP    = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_sineTrain_2Na     = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_sineTrain_Default = nan(T_isi,NP,NL,NH);
fits_time.ISIs.glm_flatTrain_3AHP    = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_flatTrain_2Na     = nan(T_isi,NP,NL,NH);
% fits_time.ISIs.glm_flatTrain_Default = nan(T_isi,NP,NL,NH);


BUFFER_d = BUFFER/downsampleRate;

dt = h*downsampleRate./1000;

refrac = 1:10;


%sts_glm(L_test,P_test,P_fit,L_fit,H,T);
%T = 4: "sine_2Na"
%T = 5: "sine_3AHP"
%T = 6: "sine_Default"

for pp_test = 1:NP
    fprintf('pp = %d / %d\n',pp_test,NP);
    TC = StimPeriods(pp_test)*1000;
    NC = floor(TT/TC);
    %tts = (1:(NC*TC))+T_range_downsampled(1)-1;
    
    xx = linspace(0,TC,NB+1);
    
    sineGainFunc2 = @(t,Period,Height)      (sin( (((t(:)-BUFFER_d)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
    sqGainFunc2   = @(t,Period,Height) round(sin( (((t(:)-BUFFER_d)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
%     sqGainFunc2   = @(t,Period,Height) round(sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
%     sineGainFunc2 = @(t,Period,Height)      (sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;

    
    
    for ll_test = 1:NL
        %%
        fprintf('  ll = %d / %d  (averaging %d cycles)\n',ll_test,NL,NC);
        
        ll_train = 4;%ll_test;
        pp_train = NP+1;
        
        g = mean(reshape(sqGainFunc2((1:(NC*TC)) + BUFFER_d,StimPeriods(pp_test),StimLevels(ll_test)),[],NC),2);
        h_Gain = zeros(NB+1,1);
        tt = 0:(TC-1);
        NB2 = nan(NB,1);
        for jj = 1:(length(xx)-1)
            h_Gain(jj,1) = sum(g(tt >= xx(jj) & tt < xx(jj+1)));
            NB2(jj) = sum(tt >= xx(jj) & tt < xx(jj+1)).*dt;
        end
        fits_time.cycleAvg_Gain(:,pp_test,ll_test) = mean(h_Gain(1:NB,:),2)./(NB2./dt);
        
        
        
%         h_glm_2Na = zeros(NB+1,NC,NH);
%         h_glm_Default = zeros(NB+1,NC,NH);
        h_glm_3AHP = zeros(NB+1,NC,NH);
        
%         sts_glm_2NA     = squeeze(sts_glm_sq_2Na(ll_test,pp_test,pp_train,ll_train,:));
        sts_glm_3AHP    = squeeze(sts_glm_sq_3AHP(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_Default = squeeze(sts_glm_sq_Default(ll_test,pp_test,pp_train,ll_train,:));
        
%         h_glm_sineTrain_2Na = zeros(NB+1,NC,NH);
%         h_glm_sineTrain_Default = zeros(NB+1,NC,NH);
        h_glm_sineTrain_3AHP = zeros(NB+1,NC,NH);
        
%         h_glm_flatTrain_2Na = zeros(NB+1,NC,NH);
%         h_glm_flatTrain_Default = zeros(NB+1,NC,NH);
        h_glm_flatTrain_3AHP = zeros(NB+1,NC,NH);
        
%         sts_glm_sineTrain_2NA     = squeeze(sts_glm_sq_sineTrain_2Na(ll_test,pp_test,pp_train,ll_train,:));
        sts_glm_sineTrain_3AHP    = squeeze(sts_glm_sq_sineTrain_3AHP(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_sineTrain_Default = squeeze(sts_glm_sq_sineTrain_Default(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_flatTrain_2NA     = squeeze(sts_glm_sq_2Na(ll_test,pp_test,pp_train,1,:));
        sts_glm_flatTrain_3AHP    = squeeze(sts_glm_sq_3AHP(ll_test,pp_test,pp_train,1,:));
%         sts_glm_flatTrain_Default = squeeze(sts_glm_sq_Default(ll_test,pp_test,pp_train,1,:));
        
%         s2Na = sts_sq_2Na{pp_test,ll_test};
        s3AHP = sts_sq_3AHP{pp_test,ll_test};
        sDefault = sts_sq_Default{pp_test,ll_test};
        
        
        
        xxs = nan(length(xx),NC);
        for cc = 1:NC
            xxs(:,cc) = xx + (cc-1)*TC + T_range_downsampled(1);
        end
        xxs = xxs(:);
%         h_2Na     = reshape(histc(s2Na,xxs),[],NC);
        h_3AHP    = reshape(histc(s3AHP,xxs),[],NC);
        h_Default = reshape(histc(sDefault,xxs),[],NC);
        
%         hspk_2Na = nan(NB,NH);
        hspk_3AHP = nan(NB,NH);
%         hspk_Default = nan(NB,NH);
%         hspkv_2Na = nan(NB,NH);
        hspkv_3AHP = nan(NB,NH);
%         hspkv_Default = nan(NB,NH);
        
        
%         hspk2_2Na = nan(NB,NH);
        hspk2_3AHP = nan(NB,NH);
%         hspk2_Default = nan(NB,NH);
%         hspk2v_2Na = nan(NB,NH);
        hspk2v_3AHP = nan(NB,NH);
%         hspk2v_Default = nan(NB,NH);
        
        
%         kstim_2Na = nan(NB,NH);
        kstim_3AHP = nan(NB,NH);
%         kstim_Default = nan(NB,NH);
%         kstimv_2Na = nan(NB,NH);
        kstimv_3AHP = nan(NB,NH);
%         kstimv_Default = nan(NB,NH);
        
        xl = round(xx(2)-xx(1));
        
%         Y_2Na_1 = zeros(length(X1),1);
%         Y_2Na_1(sts_glm_2NA{1}) = 1;
        
        Y_3AHP_1 = zeros(length(X1),1);
        Y_3AHP_1(sts_glm_3AHP{1}) = 1;
        
%         Y_Default_1 = zeros(length(X1),1);
%         Y_Default_1(sts_glm_Default{1}) = 1;
%         
        %%
        
%         fits_time.ISIs.hh_2Na(:,pp_test,ll_test) = histc(diff(s2Na),xx_isi);
        fits_time.ISIs.hh_3AHP(:,pp_test,ll_test) = histc(diff(s3AHP),xx_isi);
        fits_time.ISIs.hh_Default(:,pp_test,ll_test) = histc(diff(sDefault),xx_isi);
        for hh = 1:NH
            
%             fits_time.ISIs.glm_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_2NA{hh}),xx_isi);
            if(length(sts_glm_3AHP{hh}) > 1)
                fits_time.ISIs.glm_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_3AHP{hh}),xx_isi);
            end
%             fits_time.ISIs.glm_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_Default{hh}),xx_isi);
            
%             fits_time.ISIs.glm_sineTrain_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sineTrain_2NA{hh}),xx_isi);
            if(length(sts_glm_sineTrain_3AHP{hh}) > 1)
                fits_time.ISIs.glm_sineTrain_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sineTrain_3AHP{hh}),xx_isi);
            end
%             fits_time.ISIs.glm_sineTrain_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sineTrain_Default{hh}),xx_isi);
            
%             fits_time.ISIs.glm_flatTrain_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_2NA{hh}),xx_isi);
            if(length(sts_glm_flatTrain_3AHP{hh}) > 1)
                fits_time.ISIs.glm_flatTrain_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_3AHP{hh}),xx_isi);
            end
%             fits_time.ISIs.glm_flatTrain_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_Default{hh}),xx_isi);
        end
        
        %%
        X1_c = X1.* sqGainFunc2((1:length(X1))',StimPeriods(pp_test),StimLevels(ll_test));
        
        
        
        for hh = 1:NH
            if(isempty(sts_glm_3AHP{hh}))
                sts_glm_3AHP{hh} = [1];
            end
            if(isempty(sts_glm_sineTrain_3AHP{hh}))
                sts_glm_sineTrain_3AHP{hh} = [1];
            end
            if(isempty(sts_glm_flatTrain_3AHP{hh}))
                sts_glm_flatTrain_3AHP{hh} = [1];
            end
            
%             h_glm_2Na(:,:,hh)     = reshape(histc(sts_glm_2NA{hh},xxs),[],NC);
            h_glm_3AHP(:,:,hh)    = reshape(histc(sts_glm_3AHP{hh},xxs),[],NC);
%             h_glm_Default(:,:,hh) = reshape(histc(sts_glm_Default{hh},xxs),[],NC);
            
            
            
%             h_glm_sineTrain_2Na(:,:,hh)     = reshape(histc(sts_glm_sineTrain_2NA{hh},xxs),[],NC);
            h_glm_sineTrain_3AHP(:,:,hh)    = reshape(histc(sts_glm_sineTrain_3AHP{hh},xxs),[],NC);
%             h_glm_sineTrain_Default(:,:,hh) = reshape(histc(sts_glm_sineTrain_Default{hh},xxs),[],NC);
%             h_glm_flatTrain_2Na(:,:,hh)     = reshape(histc(sts_glm_flatTrain_2NA{hh},xxs),[],NC);
            h_glm_flatTrain_3AHP(:,:,hh)    = reshape(histc(sts_glm_flatTrain_3AHP{hh},xxs),[],NC);
%             h_glm_flatTrain_Default(:,:,hh) = reshape(histc(sts_glm_flatTrain_Default{hh},xxs),[],NC);
            
            
%             Y_c = zeros(length(X1),1);
%             Y_c(sts_glm_2NA{hh}) = 1;
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,4);
%             h_spk_c(refrac) = 0;
%             hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
%             
%             hspk_2Na(:,hh) = mean(hist_c,2);
%             hspkv_2Na(:,hh) = var(hist_c,[],2);
            
            
            Y_c = zeros(length(X1),1);
            Y_c(sts_glm_3AHP{hh}) = 1;
            h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,5);
            h_spk_c(refrac) = 0;
            hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
            hist_c = hist_c(1+(1:length(X1)));
            hist_c = reshape(hist_c(round(xxs)),[],NC);
            hist_c = hist_c(2:end,:);
            
            hspk_3AHP(:,hh) = mean(hist_c,2);
            hspkv_3AHP(:,hh) = var(hist_c,[],2);
            
            
%             Y_c = zeros(length(X1),1);
%             Y_c(sts_glm_Default{hh}) = 1;
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,6);
%             h_spk_c(refrac) = 0;
%             hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
%             
%             hspk_Default(:,hh) = mean(hist_c,2);
%             hspkv_Default(:,hh) = var(hist_c,[],2);
            
            %%
            
            
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,4);
%             h_spk_c(1:10) = 0;
%             hist_c = conv(Y_2Na_1,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
%             
%             hspk2_2Na(:,hh) = mean(hist_c,2);
%             hspk2v_2Na(:,hh) = var(hist_c,[],2);
            
            
            h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,5);
            h_spk_c(1:10) = 0;
            hist_c = conv(Y_3AHP_1,conv(h_spk_c,ones(xl,1)./xl));
            hist_c = hist_c(1+(1:length(X1)));
            hist_c = reshape(hist_c(round(xxs)),[],NC);
            hist_c = hist_c(2:end,:);
            
            hspk2_3AHP(:,hh) = mean(hist_c,2);
            hspk2v_3AHP(:,hh) = var(hist_c,[],2);
            
            
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,6);
%             h_spk_c(1:10) = 0;
%             hist_c = conv(Y_Default_1,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
%             
%             hspk2_Default(:,hh) = mean(hist_c,2);
%             hspk2v_Default(:,hh) = var(hist_c,[],2);
            
            %%
            
            
            x_stim = stimBasis*squeeze(glm_k_stim(:,pp_train,ll_train,hh,4:6)); %#ok<PFBNS>
%             st_c = conv(X1_c,x_stim(:,1));
%             st2_c = conv(st_c.^2,ones(xl,1)./xl);
%             st_c = conv(st_c,ones(xl,1)./xl);
%             st_c = st_c((1:length(X1)));
%             st_c = reshape(st_c(round(xxs)),[],NC);
%             st_c = st_c(2:end,:);
%             
%             st2_c = st2_c((1:length(X1)));
%             st2_c = reshape(st2_c(round(xxs)),[],NC);
%             st2_c = st2_c(2:end,:);
%             
%             kstim_2Na(:,hh) = mean(st_c,2);
%             kstimv_2Na(:,hh) = mean(st2_c,2) - mean(st_c,2).^2;% var(st_c,[],2);
            
            st_c = conv(X1_c,x_stim(:,2));
            st2_c = conv(st_c.^2,ones(xl,1)./xl);
            st_c = conv(st_c,ones(xl,1)./xl);
            st_c = st_c((1:length(X1)));
            st_c = reshape(st_c(round(xxs)),[],NC);
            st_c = st_c(2:end,:);
            
            st2_c = st2_c((1:length(X1)));
            st2_c = reshape(st2_c(round(xxs)),[],NC);
            st2_c = st2_c(2:end,:);
            
            kstim_3AHP(:,hh) = mean(st_c,2);
            kstimv_3AHP(:,hh) = mean(st2_c,2) - mean(st_c,2).^2;% var(st_c,[],2);
            
            
%             st_c = conv(X1_c,x_stim(:,3));
%             st2_c = conv(st_c.^2,ones(xl,1)./xl);
%             st_c = conv(st_c,ones(xl,1)./xl);
%             st_c = st_c((1:length(X1)));
%             st_c = reshape(st_c(round(xxs)),[],NC);
%             st_c = st_c(2:end,:);
%             
%             st2_c = st2_c((1:length(X1)));
%             st2_c = reshape(st2_c(round(xxs)),[],NC);
%             st2_c = st2_c(2:end,:);
%             
%             kstim_Default(:,hh) = mean(st_c,2);
%             kstim2v_Default(:,hh) = mean(st2_c,2) - mean(st_c,2).^2;% var(st_c,[],2);
        end
%         fits_time.cycleAvg_hspk_2Na(:,pp_test,ll_test,:) = hspk_2Na;
        fits_time.cycleAvg_hspk_3AHP(:,pp_test,ll_test,:) = hspk_3AHP;
%         fits_time.cycleAvg_hspk_Default(:,pp_test,ll_test,:) = hspk_Default;
        
%         fits_time.cycleAvg_hspk2_2Na(:,pp_test,ll_test,:) = hspk2_2Na;
        fits_time.cycleAvg_hspk2_3AHP(:,pp_test,ll_test,:) = hspk2_3AHP;
%         fits_time.cycleAvg_hspk2_Default(:,pp_test,ll_test,:) = hspk2_Default;
        
        
%         fits_time.cycleAvg_kstim_2Na(:,pp_test,ll_test,:) = kstim_2Na;
        fits_time.cycleAvg_kstim_3AHP(:,pp_test,ll_test,:) = kstim_3AHP;
%         fits_time.cycleAvg_kstim_Default(:,pp_test,ll_test,:) = kstim_Default;
        
%         fits_time.cycleAvg_kstimv_2Na(:,pp_test,ll_test,:) = kstimv_2Na;
        fits_time.cycleAvg_kstimv_3AHP(:,pp_test,ll_test,:) = kstimv_3AHP;
%         fits_time.cycleAvg_kstimv_Default(:,pp_test,ll_test,:) = kstimv_Default;
        
%         fits_time.cycleVar_hspk_2Na(:,pp_test,ll_test,:) = hspkv_2Na;
        fits_time.cycleVar_hspk_3AHP(:,pp_test,ll_test,:) = hspkv_3AHP;
%         fits_time.cycleVar_hspk_Default(:,pp_test,ll_test,:) = hspkv_Default;
        
        
%         fits_time.cycleAvg_2Na(:,pp_test,ll_test) = mean(h_2Na(1:NB,:),2)./(NB2);
        fits_time.cycleAvg_3AHP(:,pp_test,ll_test) = mean(h_3AHP(1:NB,:),2)./(NB2);
        fits_time.cycleAvg_Default(:,pp_test,ll_test) = mean(h_Default(1:NB,:),2)./(NB2);
        
        
%         fits_time.cycleVar_2Na(:,pp_test,ll_test) = var(h_2Na(1:NB,:)./NB2,[],2);
        fits_time.cycleVar_3AHP(:,pp_test,ll_test) = var(h_3AHP(1:NB,:)./NB2,[],2);
%         fits_time.cycleVar_Default(:,pp_test,ll_test) = var(h_Default(1:NB,:)./NB2,[],2);
        
        
%         fits_time.cycleAvg_glm_2Na(:    ,pp_test,ll_test,:) = mean(h_glm_2Na(    1:NB,:,:),2)./(NB2);
        fits_time.cycleAvg_glm_3AHP(:   ,pp_test,ll_test,:) = mean(h_glm_3AHP(   1:NB,:,:),2)./(NB2);
%         fits_time.cycleAvg_glm_Default(:,pp_test,ll_test,:) = mean(h_glm_Default(1:NB,:,:),2)./(NB2);
        
%         fits_time.cycleVar_glm_2Na(:    ,pp_test,ll_test,:) = var(h_glm_2Na(    1:NB,:,:)./NB2,[],2);
        fits_time.cycleVar_glm_3AHP(:   ,pp_test,ll_test,:) = var(h_glm_3AHP(   1:NB,:,:)./NB2,[],2);
%         fits_time.cycleVar_glm_Default(:,pp_test,ll_test,:) = var(h_glm_Default(1:NB,:,:)./NB2,[],2);
        
        
%         fits_time.cycleAvg_glm_sineTrain_2Na(:    ,pp_test,ll_test,:) = mean(h_glm_sineTrain_2Na(    1:NB,:,:),2)./(NB2);
        fits_time.cycleAvg_glm_sineTrain_3AHP(:   ,pp_test,ll_test,:) = mean(h_glm_sineTrain_3AHP(   1:NB,:,:),2)./(NB2);
%         fits_time.cycleAvg_glm_sineTrain_Default(:,pp_test,ll_test,:) = mean(h_glm_sineTrain_Default(1:NB,:,:),2)./(NB2);
        
%         fits_time.cycleVar_glm_sineTrain_2Na(:    ,pp_test,ll_test,:) = var(h_glm_sineTrain_2Na(    1:NB,:,:)./NB2,[],2);
        fits_time.cycleVar_glm_sineTrain_3AHP(:   ,pp_test,ll_test,:) = var(h_glm_sineTrain_3AHP(   1:NB,:,:)./NB2,[],2);
%         fits_time.cycleVar_glm_sineTrain_Default(:,pp_test,ll_test,:) = var(h_glm_sineTrain_Default(1:NB,:,:)./NB2,[],2);
        
%         fits_time.cycleAvg_glm_flatTrain_2Na(:    ,pp_test,ll_test,:) = mean(h_glm_flatTrain_2Na(    1:NB,:,:),2)./(NB2);
        fits_time.cycleAvg_glm_flatTrain_3AHP(:   ,pp_test,ll_test,:) = mean(h_glm_flatTrain_3AHP(   1:NB,:,:),2)./(NB2);
%         fits_time.cycleAvg_glm_flatTrain_Default(:,pp_test,ll_test,:) = mean(h_glm_flatTrain_Default(1:NB,:,:),2)./(NB2);
        
%         fits_time.cycleVar_glm_flatTrain_2Na(:    ,pp_test,ll_test,:) = var(h_glm_flatTrain_2Na(    1:NB,:,:)./NB2,[],2);
        fits_time.cycleVar_glm_flatTrain_3AHP(:   ,pp_test,ll_test,:) = var(h_glm_flatTrain_3AHP(   1:NB,:,:)./NB2,[],2);
%         fits_time.cycleVar_glm_flatTrain_Default(:,pp_test,ll_test,:) = var(h_glm_flatTrain_Default(1:NB,:,:)./NB2,[],2);
        
        
        for zz = 1:2
            tt = (1:(NB/2)) + (NB/2)*(zz-1);
            
%             [fits_time.tau_2Na(pp_test,ll_test,zz)    ,fits_time.gs_2Na(pp_test,ll_test,zz)     ,fits_time.dc_2Na(pp_test,ll_test,zz)    , fits_time.fitsAvg_2Na(tt,pp_test,ll_test)    ] = fitExp2(fits_time.cycleAvg_2Na(tt,pp_test,ll_test),StimPeriods(pp_test)/(NB));
            [fits_time.tau_3AHP(pp_test,ll_test,zz)   ,fits_time.gs_3AHP(pp_test,ll_test,zz)    ,fits_time.dc_3AHP(pp_test,ll_test,zz)   , fits_time.fitsAvg_3AHP(tt,pp_test,ll_test)   ] = fitExp2(fits_time.cycleAvg_3AHP(tt,pp_test,ll_test),StimPeriods(pp_test)/(NB));
            [fits_time.tau_Default(pp_test,ll_test,zz),fits_time.gs_Default(pp_test,ll_test,zz) ,fits_time.dc_Default(pp_test,ll_test,zz), fits_time.fitsAvg_Default(tt,pp_test,ll_test)] = fitExp2(fits_time.cycleAvg_Default(tt,pp_test,ll_test),StimPeriods(pp_test)/(NB));
            [fits_time.tau_Gain(pp_test,ll_test,zz)   ,fits_time.gs_Gain(pp_test,ll_test,zz)    ,fits_time.dc_Gain(pp_test,ll_test,zz)   , fits_time.fitsAvg_Gain(tt,pp_test,ll_test)   ] = fitExp2(fits_time.cycleAvg_Gain(tt,pp_test,ll_test),StimPeriods(pp_test)/(NB));


            for hh = 1:NH
%                 [fits_time.tau_glm_2Na(pp_test,ll_test,hh,zz)    ,fits_time.gs_glm_2Na(pp_test,ll_test,hh,zz)    ,fits_time.dc_glm_2Na(pp_test,ll_test,hh,zz)    , fits_time.fitsAvg_glm_2Na(tt,pp_test,ll_test,hh)    ] = fitExp2(fits_time.cycleAvg_glm_2Na(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
                [fits_time.tau_glm_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.gs_glm_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.dc_glm_3AHP(pp_test,ll_test,hh,zz)   , fits_time.fitsAvg_glm_3AHP(tt,pp_test,ll_test,hh)   ] = fitExp2(fits_time.cycleAvg_glm_3AHP(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
%                 [fits_time.tau_glm_Default(pp_test,ll_test,hh,zz),fits_time.gs_glm_Default(pp_test,ll_test,hh,zz),fits_time.dc_glm_Default(pp_test,ll_test,hh,zz), fits_time.fitsAvg_glm_Default(tt,pp_test,ll_test,hh)] = fitExp2(fits_time.cycleAvg_glm_Default(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));

%                 [fits_time.tau_glm_sineTrain_2Na(pp_test,ll_test,hh,zz)    ,fits_time.gs_glm_sineTrain_2Na(pp_test,ll_test,hh,zz)    ,fits_time.dc_glm_sineTrain_2Na(pp_test,ll_test,hh,zz)    , fits_time.fitsAvg_glm_sineTrain_2Na(tt,pp_test,ll_test,hh)    ] = fitExp2(fits_time.cycleAvg_glm_sineTrain_2Na(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
                [fits_time.tau_glm_sineTrain_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.gs_glm_sineTrain_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.dc_glm_sineTrain_3AHP(pp_test,ll_test,hh,zz)   , fits_time.fitsAvg_glm_sineTrain_3AHP(tt,pp_test,ll_test,hh)   ] = fitExp2(fits_time.cycleAvg_glm_sineTrain_3AHP(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
%                 [fits_time.tau_glm_sineTrain_Default(pp_test,ll_test,hh,zz),fits_time.gs_glm_sineTrain_Default(pp_test,ll_test,hh,zz),fits_time.dc_glm_sineTrain_Default(pp_test,ll_test,hh,zz), fits_time.fitsAvg_glm_sineTrain_Default(tt,pp_test,ll_test,hh)] = fitExp2(fits_time.cycleAvg_glm_sineTrain_Default(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));

%                 [fits_time.tau_glm_flatTrain_2Na(pp_test,ll_test,hh,zz)    ,fits_time.gs_glm_flatTrain_2Na(pp_test,ll_test,hh,zz)    ,fits_time.dc_glm_flatTrain_2Na(pp_test,ll_test,hh,zz)    , fits_time.fitsAvg_glm_flatTrain_2Na(tt,pp_test,ll_test,hh)    ] = fitExp2(fits_time.cycleAvg_glm_flatTrain_2Na(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
                [fits_time.tau_glm_flatTrain_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.gs_glm_flatTrain_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.dc_glm_flatTrain_3AHP(pp_test,ll_test,hh,zz)   , fits_time.fitsAvg_glm_flatTrain_3AHP(tt,pp_test,ll_test,hh)   ] = fitExp2(fits_time.cycleAvg_glm_flatTrain_3AHP(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
%                 [fits_time.tau_glm_flatTrain_Default(pp_test,ll_test,hh,zz),fits_time.gs_glm_flatTrain_Default(pp_test,ll_test,hh,zz),fits_time.dc_glm_flatTrain_Default(pp_test,ll_test,hh,zz), fits_time.fitsAvg_glm_flatTrain_Default(tt,pp_test,ll_test,hh)] = fitExp2(fits_time.cycleAvg_glm_flatTrain_Default(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));

                
%                 [fits_time.tau_hspk_2Na(pp_test,ll_test,hh,zz)    ,fits_time.gs_hspk_2Na(pp_test,ll_test,hh,zz)    ,fits_time.dc_hspk_2Na(pp_test,ll_test,hh,zz)    , fits_time.fitsAvg_hspk_2Na(tt,pp_test,ll_test,hh)    ] = fitExp2(fits_time.cycleAvg_hspk_2Na(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
                [fits_time.tau_hspk_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.gs_hspk_3AHP(pp_test,ll_test,hh,zz)   ,fits_time.dc_hspk_3AHP(pp_test,ll_test,hh,zz)   , fits_time.fitsAvg_hspk_3AHP(tt,pp_test,ll_test,hh)   ] = fitExp2(fits_time.cycleAvg_hspk_3AHP(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
%                 [fits_time.tau_hspk_Default(pp_test,ll_test,hh,zz),fits_time.gs_hspk_Default(pp_test,ll_test,hh,zz),fits_time.dc_hspk_Default(pp_test,ll_test,hh,zz), fits_time.fitsAvg_hspk_Default(tt,pp_test,ll_test,hh)] = fitExp2(fits_time.cycleAvg_hspk_Default(tt,pp_test,ll_test,hh),StimPeriods(pp_test)/(NB));
 
            end
        end
        
    end
end

save('Results/Lundstrom/Lundstrom_timeInfo.mat','-v7.3','fits_time');
    %%
figure(1);
clf;
hh = 4;

for ii = 2:4
    subplot(1,3,ii-1)
    hold on
    p1 = plot(StimPeriods,(fits_time.tau_3AHP(:,ii,1)),'k^-');
    p2 = plot(StimPeriods,(fits_time.tau_glm_3AHP(:,ii,hh,1)),'^-');
    p3 = plot(StimPeriods,(fits_time.tau_glm_sineTrain_3AHP(:,ii,hh,1)),'^-');
    p4 = plot(StimPeriods,(fits_time.tau_glm_flatTrain_3AHP(:,ii,hh,1)),'^-');
%     p3 = plot(StimPeriods,(fits_time.tau_hspk_3AHP(:,ii,hh,1)),'^-');
    
    plot(StimPeriods,(fits_time.tau_3AHP(:,ii,2)),'v-','Color',p1.Color);
    plot(StimPeriods,(fits_time.tau_glm_3AHP(:,ii,hh,2)),'v-','Color',p2.Color);
    plot(StimPeriods,(fits_time.tau_glm_sineTrain_3AHP(:,ii,hh,2)),'v-','Color',p3.Color);
    plot(StimPeriods,(fits_time.tau_glm_flatTrain_3AHP(:,ii,hh,2)),'v-','Color',p4.Color);
%     plot(StimPeriods,(fits_time.tau_hspk_3AHP(:,ii,hh,2)),'v-','Color',p3.Color);
    hold off
end


figure(2);
clf;
for ii = 2:4
    subplot(1,3,ii-1)
    p1 = loglog((StimPeriods),(fits_time.gs_3AHP(:,ii,1)),'k^-');
    hold on
    p2 = loglog((StimPeriods),(fits_time.gs_glm_3AHP(:,ii,hh,1)),'^-');
    p3 = loglog((StimPeriods),(fits_time.gs_glm_sineTrain_3AHP(:,ii,hh,1)),'^-');
    p4 = loglog((StimPeriods),(fits_time.gs_glm_flatTrain_3AHP(:,ii,hh,1)),'^-');
%     p3 = semilogx((StimPeriods),(fits_time.gs_hspk_3AHP(:,ii,hh,1)),'^-');
    
    
    loglog((StimPeriods),(fits_time.gs_3AHP(:,ii,2)),'v-','Color',p1.Color);
    loglog((StimPeriods),(fits_time.gs_glm_3AHP(:,ii,hh,2)),'v-','Color',p2.Color);
    loglog((StimPeriods),(fits_time.gs_glm_sineTrain_3AHP(:,ii,hh,2)),'v-','Color',p3.Color);
    loglog((StimPeriods),(fits_time.gs_glm_flatTrain_3AHP(:,ii,hh,2)),'v-','Color',p4.Color);
%     semilogx((StimPeriods),(fits_time.gs_hspk_3AHP(:,ii,hh,2)),'v-','Color',p3.Color);
    hold off
end