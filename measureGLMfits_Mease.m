rng(20190322);

load('Results/Mease/MeaseSims.mat','downsampleRate','sigMuRatio','StimLevels','sigMuRatio','stim_mu','xs', 'frFunctions')
load('Results/Mease/MeaseGLMs.mat','stimBasis','spkHistBases');

fr0 = squeeze(frFunctions(xs==0,:,:));

DEBUG = false;

addpath Utils/;
addpath GLMtools;
addpath kGLM;
addGLMPaths;

if(DEBUG)
    TT_start = 1*1000;
    TT_post = 1*1000;
    TT_sim = 1*1000;
    TT = TT_sim + TT_start + TT_post;
    alphas = [0 1:0.5:2];
else
    TT_start = 4*1000;
    TT_post  = 4*1000;
    TT_sim = 32*1000;
    TT = TT_sim + TT_start + TT_post;
    %alphas = [0 1:0.5:5 10];
    alphas = [0 1:1:5 10];
end


NA = length(alphas);
NH = length(spkHistBases);
poolobj = gcp('nocreate'); delete(poolobj);
parpool(NA);

smooth_sig = 5;
ff_smooth = exp(-(-500:500).^2./(2*smooth_sig^2));
ff_smooth = ff_smooth(:)./sum(ff_smooth);


NS = 1e3;

Gs = 600:100:2000;

dt = 1e-3;

T_idx = TT_start+(1:TT_sim);
TT2 = length(T_idx);

X2 = randn(TT,1);
y2 = zeros(length(X2)*downsampleRate,1);
for ii = 1:length(X2)
    y2((1:downsampleRate) + (ii-1)*downsampleRate) = X2(ii);
end

X2_c = conv2(X2,stimBasis);
X2_c = X2_c(1:TT,:);

psthSmooth_hh  = nan(TT,4,15,15);
%psth_glm = nan(TT,4,15,15,5,NH,NA);

N_Bspks = zeros(NH,1);
for H = 1:NH
    N_Bspks(H)  = size(spkHistBases{H},2);
end
glm_k_stims = nan(size(stimBasis,2),15,15,4+1,NH,NA);
glm_h_spks  = nan(max(N_Bspks),15,15,4+1,NH,NA);
glm_bs      = nan(15,15,4+1,NH,NA);
for alphaCtr = 1:NA
    alpha = alphas(alphaCtr);
    if(alpha > 0)
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    else
        saveDir = sprintf('Results/Mease/');
    end
    if(~isfolder(saveDir))
        mkdir(saveDir);
    end
    
    load(sprintf("%sMeaseGLMs.mat",saveDir),'glm_h_spk','glm_k_stim','glm_b');

    glm_h_spks(:,:,:,:,:,alphaCtr) = glm_h_spk;
    glm_k_stims(:,:,:,:,:,alphaCtr) = glm_k_stim;
    glm_bs(:,:,:,:,alphaCtr) = glm_b;
end


freqCutoff = 64;

[wcoh,wcs,pp] = wcoherence(randn(TT,1),randn(TT,1),1e3);

freqs = pp(pp <= freqCutoff);
c_test  = nan(length(freqs),4,15,15,5,NH,NA);
p_test  = nan(length(freqs),4,15,15,5,NH,NA);
p2_test = nan(length(freqs),4,15,15,5,NH,NA);
p3_test = nan(length(freqs),4,15,15,5,NH,NA);


%%
load('Results/Mease/glmFitMetrics.mat');
if(exist('Results/Mease/glmFitMetrics_part.mat','file'))
    load('Results/Mease/glmFitMetrics_part.mat')
    nn_init = nn;
    kk_init = kk;
    hh_init = hh+1;
else
    kk_init = 1;
    nn_init = 1;
    hh_init = 1;
end
%%
for nn = nn_init:15
    fprintf('nn = %d\n',nn);
    %%

    for kk = kk_init:15
        fprintf('  kk = %d\n',kk);
        kk_clock = tic;
        if(fr0(nn,kk) > 1)
            continue;
        end

        %% simulate HH model, all stim levels
        for ii = 1:4
            psthSmooth_hh(:,ii,nn,kk) = conv(psth_hh(:,ii,nn,kk),ff_smooth,'same');
        end
        psthSmooth_hh_c = psthSmooth_hh(:,:,nn,kk);
        psth_hh_c = psth_hh(:,:,nn,kk);
        startSpikes = psth_hh_c(1:TT_start,:);
        
        ps_hh  = psthSmooth_hh_c(T_idx,:);
        
        
        Fs = 1e3;
        
        cwtx_hh = wcoherence_getWaveletTransform(psth_hh(:,1),Fs);
        for ii = 2:4
            cwtx_hh(ii) = wcoherence_getWaveletTransform(psth_hh(:,ii),Fs);
        end
        
        %% for each type of GLM
        for hh = hh_init:NH
            if(isnan(glm_bs(nn,kk,1,hh,1)))
                continue;
            end
            

            
            fprintf('    hh = %d\n',hh);
            
            hh_clock = tic;
            
            spkHistBasis_c = spkHistBases{hh};

            
            
            %%
            
            glm_k_cs = squeeze(glm_k_stims(:,nn,kk,:,hh,:));
            glm_h_cs = squeeze(glm_h_spks(1:N_Bspks(hh),nn,kk,:,hh,:));
            glm_b_cs = squeeze(glm_bs(nn,kk,:,hh,:));
            
            
            
            
            parfor aa = 1:NA

                
                devNum = mod(aa,4)+1;
                gpuDevice(devNum);
                reset(gpuDevice(devNum));
                
                startSpikes_c = startSpikes;
                ps_hh_c =ps_hh;
                cwtx_hh_c = cwtx_hh;
                StimLevels_c = StimLevels;
                
                devNum = mod(aa,4);
                alpha = alphas(aa);

                %% simulate GLMs (all 5 x 10) 
                % simulate on 60s of HH data
                % use the first 5s to initialize GLM spk hist!
                glm_k_c = glm_k_cs(:,:,aa);
                glm_h_c = glm_h_cs(:,:,aa);
                glm_b_c = glm_b_cs(:,aa)';
                
                
                for tt = 1:5
                    
                    
                    %% sim glm
                
                    GLM = struct;
                    GLM.spkHistBasisVectors = spkHistBasis_c;
                    GLM.dt = dt;
                    GLM.stimBasisVectors = stimBasis;
                    if(alpha <= 0)
                        GLM.NLtype = 1;
                    elseif(alpha == 1)
                        GLM.NLtype = 2;
                    else
                        GLM.NLtype = 3;
                        GLM.NLpow = alpha;
                    end
                    GLM.k_stim = glm_k_c(:,tt);
                    GLM.h_spk  = glm_h_c(:,tt);
                    GLM.b      = glm_b_c(tt);

                    
                     %%
                    %tic;
                    stim = (X2_c*GLM.k_stim)*StimLevels_c(:)'+GLM.b;
                    
                    h_s = GLM.spkHistBasisVectors*GLM.h_spk;
                    
                    psth_glm = kcSimGLMs(h_s,GLM.dt, alpha, stim',startSpikes_c', NS, devNum);
                    %[~,psth_glm] = simGLMs_cpu(h_s,GLM.dt,alpha,NS,stim',startSpikes');
                    psth_glm = psth_glm';
                    
                    
                    %toc
                    
                    %%
%                     psthSmooth_glm = zeros(size(psth_glm));
%                     for ii = 1:4
%                         psthSmooth_glm(:,ii) = conv(psth_glm(:,ii),ff_smooth,'same');
%                     end
                    psthSmooth_glm = conv2(psth_glm,ff_smooth,'same');

                    %% compare PSTH R^2
                    
                    
                    ps_glm = psthSmooth_glm(T_idx,:);
                    
                    r2s = nan(4,1);
                    rs = nan(4,1);
                    for ii = 1:4
                        [r2s(ii),rs(ii)] = getR2(ps_hh_c(:,ii),ps_glm(:,ii));
                        
                    end
                    R2_test(:,nn,kk,tt,hh,aa) = r2s;
                    R_test(:,nn,kk,tt,hh,aa) = rs;

                    %% compare coherence
                    
                    mcohs  = nan(length(freqs),4);
                    mcss  = nan(length(freqs),4);
                    mcs2s = nan(length(freqs),4);
                    mcs3s = nan(length(freqs),4);
                    for ii = 1:4
                        Fs = 1e3;
%                         [Cxy,~] = mscohere(ps_hh(:,ii),ps_glm(:,ii),16*Fs,12*Fs,freqs(:),Fs);
%                         [Pxy,F] = cpsd(    ps_hh(:,ii),ps_glm(:,ii),16*Fs,12*Fs,freqs(:),Fs);
                        [wcoh,wcs,pp] = wcoherence_alt(cwtx_hh_c(ii),psth_glm(:,ii),Fs);
                        
                        mcohs(:,ii) = mean(real(wcoh(pp<=freqCutoff,T_idx)),2);
                        mcss(:,ii)  = mean(   (angle(wcs(pp<=freqCutoff,T_idx))),2);
                        mcs2s(:,ii) = mean(abs(angle(wcs(pp<=freqCutoff,T_idx))),2);
                        mcs3s(:,ii) = mean(abs(     (wcs(pp<=freqCutoff,T_idx))),2);
                    end
                    
                    c_test(:,:,nn,kk,tt,hh,aa) = mcohs;%Cxy(1:size(c_test,1));
                    p_test(:,:,nn,kk,tt,hh,aa) = mcss;%angle(Pxy(1:size(p_test,1)));
                    p2_test(:,:,nn,kk,tt,hh,aa) = mcs2s;
                    p3_test(:,:,nn,kk,tt,hh,aa) = mcs3s;
                    
                end
                    
                

                    
            end
            
            
            
            if(~DEBUG && (mod(hh,8) == 0 || hh == NH))
                save('Results/Mease/glmFitMetrics_part.mat','pseudo_R2_test','nn','kk','hh','alphas','X2','T_idx','c_test','p_test','p2_test','p3_test','ll_test','R_test','R2_test','freqs','psth_hh','smooth_sig','bitsPerSpike_test','-v7.3');
            end
            
            hh_time = toc(hh_clock);
            fprintf('    time for one hh set: %.2f\n',hh_time);
        end
        hh_init = 1;
        kk_time = toc(kk_clock);
        fprintf('  time for one kk set: %.2f\n',kk_time);
    end
    kk_init = 1;
end

if(~DEBUG)
    save('Results/Mease/glmFitMetrics2.mat','alphas','X2','pseudo_R2_test','T_idx','c_test','p_test','p2_test','p3_test','ll_test','R_test','R2_test','freqs','psth_hh','smooth_sig','bitsPerSpike_test','-v7.3');
end
