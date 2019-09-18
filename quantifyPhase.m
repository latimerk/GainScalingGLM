if(~exist('Results/Lundstrom/sts_glm_sine_2Na','var'))
    %load('Results/Lundstrom/LundstromGLMs.mat','sts_glm_sine_2Na','sts_glm_sine_3AHP','sts_glm_sine_Default','glm_ll','glm_h_spk','glm_k_stim','glm_b','simTypes','spkHistBases','spkHistBases_0','spkHistLengths','stimBasis');
    %load('Results/Lundstrom/LundstromGLMs_part4.mat')
    load('Results/Lundstrom/LundstromGLMs_part5.mat')
    %load('Results/Lundstrom/LundstromGLMs_part6.mat')
    %load('Results/Lundstrom/LundstromGLMs_part7.mat')
    load('Results/Lundstrom/LundstromGLMs_part8.mat')
    %load('Results/Lundstrom/LundstromGLMs_part9.mat')
    load('Results/Lundstrom/LundstromGLMs_meta.mat')
end
if(~exist('sts_sine_2Na','var'))
    load('Results/Lundstrom/LundstromSims.mat');
end

NB = 30;

NP = length(StimPeriods);
NL = length(StimLevels);
NH = length(spkHistLengths);

TT = T_range_downsampled(2)-T_range_downsampled(1)+1;

fits_phase = struct();

%fits_phase.cycleAvg_2Na     = nan(NB,NP,NL);
fits_phase.cycleAvg_3AHP    = nan(NB,NP,NL);
fits_phase.cycleAvg_Default = nan(NB,NP,NL);
fits_phase.cycleAvg_Gain    = nan(NB,NP,NL);

%fits_phase.cycleAvg_glm_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_glm_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_glm_Default = nan(NB,NP,NL,NH);

%fits_phase.cycleAvg_glm_sqTrain_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_glm_sqTrain_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_glm_sqTrain_Default = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_glm_flatTrain_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_glm_flatTrain_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_glm_flatTrain_Default = nan(NB,NP,NL,NH);


%fits_phase.cycleAvg_hspk_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_hspk_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_hspk_Default = nan(NB,NP,NL,NH);


%fits_phase.cycleAvg_kstim_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_kstim_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_kstim_Default = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_kstimv_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleAvg_kstimv_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleAvg_kstimv_Default = nan(NB,NP,NL,NH);

%fits_phase.cycleVar_hspk_2Na     = nan(NB,NP,NL,NH);
fits_phase.cycleVar_hspk_3AHP    = nan(NB,NP,NL,NH);
%fits_phase.cycleVar_hspk_Default = nan(NB,NP,NL,NH);

%fits_phase.fitsAvg_2Na     = nan(NB,NP,NL);
fits_phase.fitsAvg_3AHP    = nan(NB,NP,NL);
fits_phase.fitsAvg_Default = nan(NB,NP,NL);
fits_phase.fitsAvg_Gain    = nan(NB,NP,NL);


%fits_phase.fitsAvg_glm_2Na     = nan(NB,NP,NL,NH);
fits_phase.fitsAvg_glm_3AHP    = nan(NB,NP,NL,NH);
% fits_phase.fitsAvg_glm_Default = nan(NB,NP,NL,NH);
% fits_phase.fitsAvg_glm_sqTrain_2Na     = nan(NB,NP,NL,NH);
fits_phase.fitsAvg_glm_sqTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_phase.fitsAvg_glm_sqTrain_Default = nan(NB,NP,NL,NH);
% fits_phase.fitsAvg_glm_flatTrain_2Na     = nan(NB,NP,NL,NH);
fits_phase.fitsAvg_glm_flatTrain_3AHP    = nan(NB,NP,NL,NH);
% fits_phase.fitsAvg_glm_flatTrain_Default = nan(NB,NP,NL,NH);


% fits_phase.phase_2Na     = nan(NP,NL);
fits_phase.phase_3AHP    = nan(NP,NL);
fits_phase.phase_Default = nan(NP,NL);
fits_phase.phase_Gain    = nan(NP,NL);


% fits_phase.gs_2Na     = nan(NP,NL);
fits_phase.gs_glm_3AHP    = nan(NP,NL);
fits_phase.gs_Default = nan(NP,NL);
fits_phase.gs_Gain    = nan(NP,NL);

% fits_phase.dc_2Na     = nan(NP,NL);
fits_phase.dc_3AHP    = nan(NP,NL);
fits_phase.dc_Default = nan(NP,NL);
fits_phase.dc_Gain    = nan(NP,NL);


 
% fits_phase.phase_glm_2Na     = nan(NP,NL,NH);
fits_phase.phase_glm_3AHP    = nan(NP,NL,NH);
% fits_phase.phase_glm_Default = nan(NP,NL,NH);

% fits_phase.gs_glm_2Na     = nan(NP,NL,NH);
fits_phase.gs_glm_3AHP    = nan(NP,NL,NH);
% fits_phase.gs_glm_Default = nan(NP,NL,NH);

% fits_phase.dc_glm_2Na     = nan(NP,NL,NH);
fits_phase.dc_glm_3AHP    = nan(NP,NL,NH);
% fits_phase.dc_glm_Default = nan(NP,NL,NH);

% fits_phase.phase_glm_sqTrain_2Na     = nan(NP,NL,NH);
fits_phase.phase_glm_sqTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.phase_glm_sqTrain_Default = nan(NP,NL,NH);
% fits_phase.phase_glm_flatTrain_2Na     = nan(NP,NL,NH);
fits_phase.phase_glm_flatTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.phase_glm_flatTrain_Default = nan(NP,NL,NH);

% fits_phase.gs_glm_sqTrain_2Na     = nan(NP,NL,NH);
fits_phase.gs_glm_sqTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.gs_glm_sqTrain_Default = nan(NP,NL,NH);
% fits_phase.gs_glm_flatTrain_2Na     = nan(NP,NL,NH);
fits_phase.gs_glm_flatTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.gs_glm_flatTrain_Default = nan(NP,NL,NH);

% fits_phase.dc_glm_sqTrain_2Na     = nan(NP,NL,NH);
fits_phase.dc_glm_sqTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.dc_glm_sqTrain_Default = nan(NP,NL,NH);
% fits_phase.dc_glm_flatTrain_2Na     = nan(NP,NL,NH);
fits_phase.dc_glm_flatTrain_3AHP    = nan(NP,NL,NH);
% fits_phase.dc_glm_flatTrain_Default = nan(NP,NL,NH);



% fits_phase.phase_hspk_2Na     = nan(NP,NL,NH);
fits_phase.phase_hspk_3AHP    = nan(NP,NL,NH);
% fits_phase.phase_hspk_Default = nan(NP,NL,NH);

% fits_phase.gs_hspk_2Na     = nan(NP,NL,NH);
fits_phase.gs_hspk_3AHP    = nan(NP,NL,NH);
% fits_phase.gs_hspk_Default = nan(NP,NL,NH);

% fits_phase.dc_hspk_2Na     = nan(NP,NL,NH);
fits_phase.dc_hspk_3AHP    = nan(NP,NL,NH);
% fits_phase.dc_hspk_Default = nan(NP,NL,NH);

BUFFER_d = BUFFER/downsampleRate;

dt = h*downsampleRate./1000;


xx_isi = 0:1:1000;
T_isi = length(xx_isi);
fits_phase.ISIs.xx = xx_isi;
fits_phase.ISIs.hh_3AHP    = nan(T_isi,NP,NL);
% fits_phase.ISIs.hh_2Na     = nan(T_isi,NP,NL);
fits_phase.ISIs.hh_Default = nan(T_isi,NP,NL);
fits_phase.ISIs.glm_3AHP   = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_2Na    = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_Default = nan(T_isi,NP,NL,NH);
fits_phase.ISIs.glm_sqTrain_3AHP    = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_sqTrain_2Na     = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_sqTrain_Default = nan(T_isi,NP,NL,NH);
fits_phase.ISIs.glm_flatTrain_3AHP    = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_flatTrain_2Na     = nan(T_isi,NP,NL,NH);
% fits_phase.ISIs.glm_flatTrain_Default = nan(T_isi,NP,NL,NH);



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

    
    pp_train = NP+1;
    
    for ll_test = 1:NL
        ll_train = 4;%ll_test;
        fprintf('  ll = %d / %d\n',ll_test,NL);
        
        g = mean(reshape(sineGainFunc2((1:(NC*TC)) + BUFFER_d,StimPeriods(pp_test),StimLevels(ll_test)),[],NC),2);
        h_Gain = zeros(NB+1,1);
        tt = 0:(TC-1);
        NB2 = nan(NB,1);
        for jj = 1:(length(xx)-1)
            h_Gain(jj,1) = sum(g(tt >= xx(jj) & tt < xx(jj+1)));
            NB2(jj) = sum(tt >= xx(jj) & tt < xx(jj+1)).*dt;
        end
        fits_phase.cycleAvg_Gain(:,pp_test,ll_test) = mean(h_Gain(1:NB,:),2)./(NB2./dt);
        
        
%         h_glm_2Na = zeros(NB+1,NC,NH);
        h_glm_Default = zeros(NB+1,NC,NH);
        h_glm_3AHP = zeros(NB+1,NC,NH);
        
        
%         h_glm_sqTrain_2Na = zeros(NB+1,NC,NH);
%         h_glm_sqTrain_Default = zeros(NB+1,NC,NH);
        h_glm_sqTrain_3AHP = zeros(NB+1,NC,NH);
%         h_glm_flatTrain_2Na = zeros(NB+1,NC,NH);
%         h_glm_flatTrain_Default = zeros(NB+1,NC,NH);
        h_glm_flatTrain_3AHP = zeros(NB+1,NC,NH);
        
%         sts_glm_2NA     = squeeze(sts_glm_sine_2Na(ll_test,pp_test,pp_train,ll_train,:));
        sts_glm_3AHP    = squeeze(sts_glm_sine_3AHP(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_Default = squeeze(sts_glm_sine_Default(ll_test,pp_test,pp_train,ll_train,:));
        
%         sts_glm_sqTrain_2NA     = squeeze(sts_glm_sine_sqTrain_2Na(ll_test,pp_test,pp_train,ll_train,:));
        sts_glm_sqTrain_3AHP    = squeeze(sts_glm_sine_sqTrain_3AHP(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_sqTrain_Default = squeeze(sts_glm_sine_sqTrain_Default(ll_test,pp_test,pp_train,ll_train,:));
%         sts_glm_flatTrain_2NA     = squeeze(sts_glm_sine_2Na(ll_test,pp_test,pp_train,1,:));
        sts_glm_flatTrain_3AHP    = squeeze(sts_glm_sine_3AHP(ll_test,pp_test,pp_train,1,:));
%         sts_glm_flatTrain_Default = squeeze(sts_glm_sine_Default(ll_test,pp_test,pp_train,1,:));
        
%         s2Na = sts_sine_2Na{pp_test,ll_test};
        s3AHP = sts_sine_3AHP{pp_test,ll_test};
        sDefault = sts_sine_Default{pp_test,ll_test};
        
%         fprintf('    cc = %d / %d\n',cc,NC);
%         h_2Na = zeros(NB+1,NC);
%         h_Default = zeros(NB+1,NC);
%         h_3AHP = zeros(NB+1,NC);
%         
%         for cc = 1:NC
%             xxs = xx + (cc-1)*TC + T_range_downsampled(1);
%             h_2Na(:,cc)     = histc(s2Na,xxs);
%             h_3AHP(:,cc)    = histc(s3AHP,xxs);
%             h_Default(:,cc) = histc(sDefault,xxs);
%             
%             for hh = 1:NH
%                 h_glm_2Na(:,cc,hh)     = histc(sts_glm_2NA{hh},xxs);
%                 h_glm_3AHP(:,cc,hh)    = histc(sts_glm_3AHP{hh},xxs);
%                 h_glm_Default(:,cc,hh) = histc(sts_glm_Default{hh},xxs);
%             end
%         end

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
        
        
%         kstim_2Na = nan(NB,NH);
        kstim_3AHP = nan(NB,NH);
%         kstim_Default = nan(NB,NH);
%         kstimv_2Na = nan(NB,NH);
        kstimv_3AHP = nan(NB,NH);
%         kstimv_Default = nan(NB,NH);
        
        
        X1_c = X1.* sineGainFunc2((1:length(X1))',StimPeriods(pp_test),StimLevels(ll_test));
        
        xl = round(xx(2)-xx(1));
        %%
        
%         fits_phase.ISIs.hh_2Na(:,pp_test,ll_test) = histc(diff(s2Na),xx_isi);
        fits_phase.ISIs.hh_3AHP(:,pp_test,ll_test) = histc(diff(s3AHP),xx_isi);
        fits_phase.ISIs.hh_Default(:,pp_test,ll_test) = histc(diff(sDefault),xx_isi);
        for hh = 1:NH
            
%             fits_phase.ISIs.glm_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_2NA{hh}),xx_isi);
            fits_phase.ISIs.glm_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_3AHP{hh}),xx_isi);
%             fits_phase.ISIs.glm_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_Default{hh}),xx_isi);
            
%             fits_phase.ISIs.glm_sqTrain_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sqTrain_2NA{hh}),xx_isi);
            fits_phase.ISIs.glm_sqTrain_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sqTrain_3AHP{hh}),xx_isi);
%             fits_phase.ISIs.glm_sqTrain_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_sqTrain_Default{hh}),xx_isi);
            
%             fits_phase.ISIs.glm_flatTrain_2Na(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_2NA{hh}),xx_isi);
            fits_phase.ISIs.glm_flatTrain_3AHP(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_3AHP{hh}),xx_isi);
%             fits_phase.ISIs.glm_flatTrain_Default(:,pp_test,ll_test,hh) = histc(diff(sts_glm_flatTrain_Default{hh}),xx_isi);
        end
        %%
        for hh = 1:NH
%             h_glm_2Na(:,:,hh)     = reshape(histc(sts_glm_2NA{hh},xxs),[],NC);
            h_glm_3AHP(:,:,hh)    = reshape(histc(sts_glm_3AHP{hh},xxs),[],NC);
%             h_glm_Default(:,:,hh) = reshape(histc(sts_glm_Default{hh},xxs),[],NC);
            
%             h_glm_sqTrain_2Na(:,:,hh)     = reshape(histc(sts_glm_sqTrain_2NA{hh},xxs),[],NC);
            h_glm_sqTrain_3AHP(:,:,hh)    = reshape(histc(sts_glm_sqTrain_3AHP{hh},xxs),[],NC);
%             h_glm_sqTrain_Default(:,:,hh) = reshape(histc(sts_glm_sqTrain_Default{hh},xxs),[],NC);
%             h_glm_flatTrain_2Na(:,:,hh)     = reshape(histc(sts_glm_flatTrain_2NA{hh},xxs),[],NC);
            h_glm_flatTrain_3AHP(:,:,hh)    = reshape(histc(sts_glm_flatTrain_3AHP{hh},xxs),[],NC);
%             h_glm_flatTrain_Default(:,:,hh) = reshape(histc(sts_glm_flatTrain_Default{hh},xxs),[],NC);
            
            
%             Y_c = zeros(length(X1),1);
%             Y_c(sts_glm_2NA{hh}) = 1;
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,1);
%             h_spk_c(1:10) = 0;
%             hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
            
%             hspk_2Na(:,hh) = mean(hist_c,2);
%             hspkv_2Na(:,hh) = var(hist_c,[],2);
            
            
            Y_c = zeros(length(X1),1);
            Y_c(sts_glm_3AHP{hh}) = 1;
            h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,2);
            h_spk_c(1:10) = 0;
            hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
            hist_c = hist_c(1+(1:length(X1)));
            hist_c = reshape(hist_c(round(xxs)),[],NC);
            hist_c = hist_c(2:end,:);
            
            hspk_3AHP(:,hh) = mean(hist_c,2);
            hspkv_3AHP(:,hh) = var(hist_c,[],2);
            
            
%             Y_c = zeros(length(X1),1);
%             Y_c(sts_glm_Default{hh}) = 1;
%             h_spk_c = spkHistBases{hh}*glm_h_spk(1:size(spkHistBases{hh},2),pp_train,ll_train,hh,3);
%             h_spk_c(1:10) = 0;
%             hist_c = conv(Y_c,conv(h_spk_c,ones(xl,1)./xl));
%             hist_c = hist_c(1+(1:length(X1)));
%             hist_c = reshape(hist_c(round(xxs)),[],NC);
%             hist_c = hist_c(2:end,:);
%             
%             hspk_Default(:,hh) = mean(hist_c,2);
%             hspkv_Default(:,hh) = var(hist_c,[],2);
            
             %%
            
            
            x_stim = stimBasis*squeeze(glm_k_stim(:,pp_train,ll_train,hh,1:3)); %#ok<PFBNS>
            
            
            
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
            
 
        end
%         fits_phase.cycleAvg_hspk_2Na(:,pp_test,ll_test,:) = hspk_2Na;
        fits_phase.cycleAvg_hspk_3AHP(:,pp_test,ll_test,:) = hspk_3AHP;
%         fits_phase.cycleAvg_hspk_Default(:,pp_test,ll_test,:) = hspk_Default;
%         fits_phase.cycleVar_hspk_2Na(:,pp_test,ll_test,:) = hspkv_2Na;
        fits_phase.cycleVar_hspk_3AHP(:,pp_test,ll_test,:) = hspkv_3AHP;
%         fits_phase.cycleVar_hspk_Default(:,pp_test,ll_test,:) = hspkv_Default;
        
        
%         fits_phase.cycleAvg_2Na(:,pp_test,ll_test) = mean(h_2Na(1:NB,:),2)./(NB2);
        fits_phase.cycleAvg_3AHP(:,pp_test,ll_test) = mean(h_3AHP(1:NB,:),2)./(NB2);
        fits_phase.cycleAvg_Default(:,pp_test,ll_test) = mean(h_Default(1:NB,:),2)./(NB2);
        
        
        
%         fits_phase.cycleAvg_kstim_2Na(:,pp_test,ll_test,:) = kstim_2Na;
        fits_phase.cycleAvg_kstim_3AHP(:,pp_test,ll_test,:) = kstim_3AHP;
%         fits_phase.cycleAvg_kstim_Default(:,pp_test,ll_test,:) = kstim_Default;
        
%         fits_phase.cycleAvg_kstimv_2Na(:,pp_test,ll_test,:) = kstimv_2Na;
        fits_phase.cycleAvg_kstimv_3AHP(:,pp_test,ll_test,:) = kstimv_3AHP;
%         fits_phase.cycleAvg_kstimv_Default(:,pp_test,ll_test,:) = kstimv_Default;
        
        
%         fits_phase.cycleVar_2Na(:,pp_test,ll_test) = var(h_2Na(1:NB,:)./NB2,[],2);
        fits_phase.cycleVar_3AHP(:,pp_test,ll_test) = var(h_3AHP(1:NB,:)./NB2,[],2);
        fits_phase.cycleVar_Default(:,pp_test,ll_test) = var(h_Default(1:NB,:)./NB2,[],2);
        
        
%         fits_phase.cycleAvg_glm_2Na(:,pp_test,ll_test,:)  = mean(h_glm_2Na(1:NB,:,:),2)./(NB2);
        fits_phase.cycleAvg_glm_3AHP(:,pp_test,ll_test,:) = mean(h_glm_3AHP(1:NB,:,:),2)./(NB2);
%         fits_phase.cycleAvg_glm_Default(:,pp_test,ll_test,:) = mean(h_glm_Default(1:NB,:,:),2)./(NB2);
        
%         fits_phase.cycleVar_glm_2Na(:,pp_test,ll_test,:)  = var(h_glm_2Na(1:NB,:,:)./NB2,[],2);
        fits_phase.cycleVar_glm_3AHP(:,pp_test,ll_test,:) = var(h_glm_3AHP(1:NB,:,:)./NB2,[],2);
%         fits_phase.cycleVar_glm_Default(:,pp_test,ll_test,:) = var(h_glm_Default(1:NB,:,:)./NB2,[],2);
        
        
        
%         fits_phase.cycleAvg_glm_sqTrain_2Na(:,pp_test,ll_test,:)  = mean(h_glm_sqTrain_2Na(1:NB,:,:),2)./(NB2);
        fits_phase.cycleAvg_glm_sqTrain_3AHP(:,pp_test,ll_test,:) = mean(h_glm_sqTrain_3AHP(1:NB,:,:),2)./(NB2);
%         fits_phase.cycleAvg_glm_sqTrain_Default(:,pp_test,ll_test,:) = mean(h_glm_sqTrain_Default(1:NB,:,:),2)./(NB2);
%         fits_phase.cycleAvg_glm_flatTrain_2Na(:,pp_test,ll_test,:)  = mean(h_glm_flatTrain_2Na(1:NB,:,:),2)./(NB2);
        fits_phase.cycleAvg_glm_flatTrain_3AHP(:,pp_test,ll_test,:) = mean(h_glm_flatTrain_3AHP(1:NB,:,:),2)./(NB2);
%         fits_phase.cycleAvg_glm_flatTrain_Default(:,pp_test,ll_test,:) = mean(h_glm_flatTrain_Default(1:NB,:,:),2)./(NB2);
        
%         fits_phase.cycleVar_glm_sqTrain_2Na(:,pp_test,ll_test,:)  = var(h_glm_sqTrain_2Na(1:NB,:,:)./NB2,[],2);
        fits_phase.cycleVar_glm_sqTrain_3AHP(:,pp_test,ll_test,:) = var(h_glm_sqTrain_3AHP(1:NB,:,:)./NB2,[],2);
%         fits_phase.cycleVar_glm_sqTrain_Default(:,pp_test,ll_test,:) = var(h_glm_sqTrain_Default(1:NB,:,:)./NB2,[],2);
%         fits_phase.cycleVar_glm_flatTrain_2Na(:,pp_test,ll_test,:)  = var(h_glm_flatTrain_2Na(1:NB,:,:)./NB2,[],2);
        fits_phase.cycleVar_glm_flatTrain_3AHP(:,pp_test,ll_test,:) = var(h_glm_flatTrain_3AHP(1:NB,:,:)./NB2,[],2);
%         fits_phase.cycleVar_glm_flatTrain_Default(:,pp_test,ll_test,:) = var(h_glm_flatTrain_Default(1:NB,:,:)./NB2,[],2);
        
%         [fits_phase.phase_2Na(pp_test,ll_test)    ,fits_phase.gs_2Na(pp_test,ll_test)    ,fits_phase.dc_2Na(pp_test,ll_test)     , fits_phase.fitsAvg_2Na(:,pp_test,ll_test)    ] = fitSine(fits_phase.cycleAvg_2Na(:,pp_test,ll_test));
        [fits_phase.phase_3AHP(pp_test,ll_test)   ,fits_phase.gs_3AHP(pp_test,ll_test)   ,fits_phase.dc_3AHP(pp_test,ll_test)    , fits_phase.fitsAvg_3AHP(:,pp_test,ll_test)   ] = fitSine(fits_phase.cycleAvg_3AHP(:,pp_test,ll_test));
        [fits_phase.phase_Default(pp_test,ll_test),fits_phase.gs_Default(pp_test,ll_test),fits_phase.dc_Default(pp_test,ll_test) , fits_phase.fitsAvg_Default(:,pp_test,ll_test)] = fitSine(fits_phase.cycleAvg_Default(:,pp_test,ll_test));
        [fits_phase.phase_Gain(pp_test,ll_test)   ,fits_phase.gs_Gain(pp_test,ll_test)   ,fits_phase.dc_Gain(pp_test,ll_test)    , fits_phase.fitsAvg_Gain(:,pp_test,ll_test)   ] = fitSine(fits_phase.cycleAvg_Gain(:,pp_test,ll_test));
        
        for hh = 1:NH
%             [fits_phase.phase_glm_2Na(pp_test,ll_test,hh)    ,fits_phase.gs_glm_2Na(pp_test,ll_test,hh)    ,fits_phase.dc_glm_2Na(pp_test,ll_test,hh)    , fits_phase.fitsAvg_glm_2Na(:,pp_test,ll_test,hh)    ] = fitSine(fits_phase.cycleAvg_glm_2Na(:,pp_test,ll_test,hh));
            [fits_phase.phase_glm_3AHP(pp_test,ll_test,hh)   ,fits_phase.gs_glm_3AHP(pp_test,ll_test,hh)   ,fits_phase.dc_glm_3AHP(pp_test,ll_test,hh)   , fits_phase.fitsAvg_glm_3AHP(:,pp_test,ll_test,hh)   ] = fitSine(fits_phase.cycleAvg_glm_3AHP(:,pp_test,ll_test,hh));
%             [fits_phase.phase_glm_Default(pp_test,ll_test,hh),fits_phase.gs_glm_Default(pp_test,ll_test,hh),fits_phase.dc_glm_Default(pp_test,ll_test,hh), fits_phase.fitsAvg_glm_Default(:,pp_test,ll_test,hh)] = fitSine(fits_phase.cycleAvg_glm_Default(:,pp_test,ll_test,hh));
        
%             [fits_phase.phase_glm_sqTrain_2Na(pp_test,ll_test,hh)    ,fits_phase.gs_glm_sqTrain_2Na(pp_test,ll_test,hh)    ,fits_phase.dc_glm_sqTrain_2Na(pp_test,ll_test,hh)    , fits_phase.fitsAvg_glm_sqTrain_2Na(:,pp_test,ll_test,hh)    ] = fitSine(fits_phase.cycleAvg_glm_sqTrain_2Na(:,pp_test,ll_test,hh));
            [fits_phase.phase_glm_sqTrain_3AHP(pp_test,ll_test,hh)   ,fits_phase.gs_glm_sqTrain_3AHP(pp_test,ll_test,hh)   ,fits_phase.dc_glm_sqTrain_3AHP(pp_test,ll_test,hh)   , fits_phase.fitsAvg_glm_sqTrain_3AHP(:,pp_test,ll_test,hh)   ] = fitSine(fits_phase.cycleAvg_glm_sqTrain_3AHP(:,pp_test,ll_test,hh));
%             [fits_phase.phase_glm_sqTrain_Default(pp_test,ll_test,hh),fits_phase.gs_glm_sqTrain_Default(pp_test,ll_test,hh),fits_phase.dc_glm_sqTrain_Default(pp_test,ll_test,hh), fits_phase.fitsAvg_glm_sqTrain_Default(:,pp_test,ll_test,hh)] = fitSine(fits_phase.cycleAvg_glm_sqTrain_Default(:,pp_test,ll_test,hh));
%             [fits_phase.phase_glm_flatTrain_2Na(pp_test,ll_test,hh)    ,fits_phase.gs_glm_flatTrain_2Na(pp_test,ll_test,hh)    ,fits_phase.dc_glm_flatTrain_2Na(pp_test,ll_test,hh)    , fits_phase.fitsAvg_glm_flatTrain_2Na(:,pp_test,ll_test,hh)    ] = fitSine(fits_phase.cycleAvg_glm_flatTrain_2Na(:,pp_test,ll_test,hh));
            [fits_phase.phase_glm_flatTrain_3AHP(pp_test,ll_test,hh)   ,fits_phase.gs_glm_flatTrain_3AHP(pp_test,ll_test,hh)   ,fits_phase.dc_glm_flatTrain_3AHP(pp_test,ll_test,hh)   , fits_phase.fitsAvg_glm_flatTrain_3AHP(:,pp_test,ll_test,hh)   ] = fitSine(fits_phase.cycleAvg_glm_flatTrain_3AHP(:,pp_test,ll_test,hh));
%             [fits_phase.phase_glm_flatTrain_Default(pp_test,ll_test,hh),fits_phase.gs_glm_flatTrain_Default(pp_test,ll_test,hh),fits_phase.dc_glm_flatTrain_Default(pp_test,ll_test,hh), fits_phase.fitsAvg_glm_flatTrain_Default(:,pp_test,ll_test,hh)] = fitSine(fits_phase.cycleAvg_glm_flatTrain_Default(:,pp_test,ll_test,hh));
        
            
%             [fits_phase.phase_hspk_2Na(pp_test,ll_test,hh)    ,fits_phase.gs_hspk_2Na(pp_test,ll_test,hh)    ,fits_phase.dc_hspk_2Na(pp_test,ll_test,hh)    , fits_phase.fitsAvg_hspk_2Na(:,pp_test,ll_test,hh)    ] = fitSine(fits_phase.cycleAvg_hspk_2Na(:,pp_test,ll_test,hh));
            [fits_phase.phase_hspk_3AHP(pp_test,ll_test,hh)   ,fits_phase.gs_hspk_3AHP(pp_test,ll_test,hh)   ,fits_phase.dc_hspk_3AHP(pp_test,ll_test,hh)   , fits_phase.fitsAvg_hspk_3AHP(:,pp_test,ll_test,hh)   ] = fitSine(fits_phase.cycleAvg_hspk_3AHP(:,pp_test,ll_test,hh));
%             [fits_phase.phase_hspk_Default(pp_test,ll_test,hh),fits_phase.gs_hspk_Default(pp_test,ll_test,hh),fits_phase.dc_hspk_Default(pp_test,ll_test,hh), fits_phase.fitsAvg_hspk_Default(:,pp_test,ll_test,hh)] = fitSine(fits_phase.cycleAvg_hspk_Default(:,pp_test,ll_test,hh));
        
        end
    end
end

save('Results/Lundstrom/Lundstrom_phaseInfo.mat','-v7.3','fits_phase');
    %%
hh_test = 4;
figure(1);
clf;
for ii = 2:4
    subplot(2,2,ii)
    hold on
%     plot(StimPeriods,rad2deg(fits_phase.phase_2Na(:,ii)-fits_phase.phase_Gain(:,ii)),'bv--')
%     plot(StimPeriods,rad2deg(fits_phase.phase_glm_2Na(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'bo-')
    plot(StimPeriods,rad2deg(fits_phase.phase_3AHP(:,ii)-fits_phase.phase_Gain(:,ii)),'ko-')
    plot(StimPeriods,rad2deg(fits_phase.phase_glm_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-')
    plot(StimPeriods,rad2deg(fits_phase.phase_glm_sqTrain_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-')
    plot(StimPeriods,rad2deg(fits_phase.phase_glm_flatTrain_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-')
%     plot(StimPeriods,rad2deg(fits_phase.phase_hspk_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o:')
%     plot(StimPeriods,rad2deg(fits_phase.phase_Default(:,ii)-fits_phase.phase_Gain(:,ii)),'gv--')
%     plot(StimPeriods,rad2deg(fits_phase.phase_glm_Default(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'go-')
    hold off
end


figure(2);
clf;
for ii = 2:4
    subplot(2,2,ii)
%     semilogx((StimPeriods),(fits_phase.gs_2Na(:,ii)./fits_phase.gs_Gain(:,ii)),'bv--')
%     hold on
%     semilogx((StimPeriods),(fits_phase.gs_glm_2Na(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'bo-')
    
    loglog((StimPeriods),(fits_phase.gs_3AHP(:,ii)./fits_phase.gs_Gain(:,ii)),'ko--')
    hold on
    loglog((StimPeriods),(fits_phase.gs_glm_3AHP(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'o-')
    loglog((StimPeriods),(fits_phase.gs_glm_sqTrain_3AHP(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'o-')
    loglog((StimPeriods),(fits_phase.gs_glm_flatTrain_3AHP(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'o-')
%     semilogx((StimPeriods),(fits_phase.gs_hspk_3AHP(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'ro:')
%     semilogx((StimPeriods),(fits_phase.gs_Default(:,ii)./fits_phase.gs_Gain(:,ii)),'gv--')
%     semilogx((StimPeriods),(fits_phase.gs_glm_Default(:,ii,hh_test)./fits_phase.gs_Gain(:,ii)),'go-')
    hold off
end
