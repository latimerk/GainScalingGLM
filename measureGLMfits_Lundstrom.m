
load('Results/Lundstrom/LundstromGLMs_meta.mat')
load('Results/Lundstrom/LundstromSims.mat', 'sigMuRatio','stim_mu_3AHP','StimLevels','downsampleRate')

addpath Utils/;
addpath GLMtools;
addpath kGLM;
addGLMPaths;

spkHistLengths(1) = 10;
rng(03262019);

%% setup glms
NH = length(spkHistBases);
NL = length(StimLevels);


T_train = 8;
L_train = 4;

glm_ks = squeeze(glm_k_stim(:,T_train,L_train,:,[5 2])); %order is sine, square, flat
glm_ks = cat(3,glm_ks,squeeze(glm_k_stim(:,T_train,1,:,5)));

glm_hs = squeeze(glm_h_spk(:,T_train,L_train,:,[5 2]));
glm_hs = cat(3,glm_hs,squeeze(glm_h_spk(:,T_train,1,:,5)));

glm_bs = squeeze(glm_b(T_train,L_train,:,[5 2]));
glm_bs = cat(2,glm_bs,squeeze(glm_b(T_train,1,:,5)));

NG = size(glm_ks,3);

%% setup stimulus


DEBUG = false;
if(DEBUG)
    TT_start = 1*1000;
    TT_post = 1*1000;
    TT_sim = 1*1000;
    TT = TT_post + TT_sim + TT_start;
else
    TT_start = 64*1000;
    TT_post = 4*1000;
    TT_sim = 32*1000;
    TT = TT_post + TT_sim + TT_start;
end
T_idx = TT_start+(1:TT_sim);

dt = 1e-3;
TT2 = length(T_idx);

X2 = randn(TT,1);

NC = 2;
conds = cell(NC,2);

h = 1/downsampleRate;
BUFFER = TT_start*downsampleRate;
Period_const = 4;

conds{1,1} = @(t)      (sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period_const).*0.5+0.5).*(StimLevels(:)'-1) + 1;
conds{1,2} = @(t,Height)      (sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period_const).*0.5+0.5).*(Height(:)'-1) + 1;

conds{2,1} = @(t) round(sin( (((t(:)-BUFFER)*h)/1000.0 *2*pi)./Period_const).*0.5+0.5).*(StimLevels(:)'-1) + 1;
conds{2,2} = @(t,Height) round(sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period_const).*0.5+0.5).*(Height(:)'-1) + 1;

y2 = zeros(length(X2)*downsampleRate,1);
for ii = 1:length(X2)
    y2((1:downsampleRate) + (ii-1)*downsampleRate,:) = X2(ii,:);
end

%% setup storage

smooth_sig = 5;
ff_smooth = exp(-(-500:500).^2./(2*smooth_sig^2));
ff_smooth = ff_smooth(:)./sum(ff_smooth);

NS = 1e3;

psth_hh        = nan(TT,NL,NC);
psthSmooth_hh  = nan(TT,NL,NC);
psth_glm       = nan(TT,NL,NC,NH,NG);

N_Bspks = zeros(NH,1);
for H = 1:NH
    N_Bspks(H)  = size(spkHistBases{H},2);
end

freqCutoff = 64;

[wcoh,wcs,pp] = wcoherence(randn(TT,1),randn(TT,1),1e3);

freqs = pp(pp <= freqCutoff);
c_test  = nan(length(freqs),NL,NC,NH,NG);
p_test  = nan(length(freqs),NL,NC,NH,NG);
p2_test = nan(length(freqs),NL,NC,NH,NG);
p3_test = nan(length(freqs),NL,NC,NH,NG);

ll_test           = nan(NL,NC,NH,NG);
bitsPerSpike_test = nan(NL,NC,NH,NG);
R_test            = nan(NL,NC,NH,NG);
R2_test           = nan(NL,NC,NH,NG);
pseudo_R2_test           = nan(NL,NC,NH,NG);

devNum = 0;

% load('Results/Lundstrom/glmFitMetrics_hhSims.mat');
figure(1);
clf;
drawnow;

%% for each condition (StimLevel and sine/square)
for cc = 1:NC
    
    
    %% simulate HH model 
    sts = simLundstrom_3AHP(y2,h,stim_mu_3AHP*sigMuRatio,stim_mu_3AHP,conds{cc,1});
    for ll = 1:NL
        st = ceil(sts{ll}./downsampleRate);
        psth_hh(: ,ll,cc) = 0;
        psth_hh(st,ll,cc) = 1;
        
        psthSmooth_hh(: ,ll,cc) = conv2(psth_hh(: ,ll,cc),ff_smooth,'same');
    end

    for ll = 1:NL %StimLevels (gains)
        fprintf('  ll = %d / %d  (cc = %d)\n',ll,NL,cc);
        %% get glm stimulus
        gain = conds{cc,2}(1:TT,StimLevels(ll));
        X3 = X2.*gain;
        X_c = conv2(X3,stimBasis);
        X_c = X_c(1:TT,:);
        
        for hh = 1:NH
            fprintf('    hh = %d / %d  (cc = %d)\n',hh,NH,cc);
            spkHistBasis_c = spkHistBases{hh};
            HT = size(spkHistBasis_c,1);
            NB = size(spkHistBasis_c,2);
            
            %% prep spike history
            Y_hist = conv2(psth_hh(:,ll,cc),spkHistBasis_c);
            Y_hist = [zeros(1,NB);Y_hist(1:(TT-1),:)];
            
            
            XX_glm = [X_c(T_idx,:) Y_hist(T_idx,:) ones(length(T_idx),1)];
            Y_c    = psth_hh(T_idx,ll,cc);
            Xy = XX_glm'*Y_c;
            normalizer = (log(dt)*sum(Y_c,1)-sum(gammaln(Y_c+1),1));
            
            %%
            for gg = 1:NG
                if(~isnan(p3_test(1,ll,cc,hh,gg)))
                    continue;
                end
                
                fprintf('      gg = %d / %d  (cc = %d)\n',gg,NG,cc);
                %% get current glm params
                glm_k = glm_ks(:,hh,gg);
                glm_h = glm_hs(1:NB,hh,gg);
                glm_b = glm_bs(hh,gg);
                glm_h_filt = spkHistBasis_c*glm_h;
                
                %% get GLM bits per spike 
                
                nllFunc = @(w) glmNll_PoissonExp(w,XX_glm,Y_c,dt,Xy,normalizer);
                %nllFunc = @(w) glmNll_BerExp(w,XX_glm,Y_c,dt);
                
                ll_model = -nllFunc([glm_k;glm_h;glm_b]);
                ll_test(ll,cc,hh,gg) =  ll_model;

                n_sp = sum(Y_c);
                T = size(Y_c,1);

                ll_null = -n_sp + n_sp*(log(n_sp) - log(T));
%                 p = n_sp/T;
%                 ll_null = (T-n_sp)*log(1-p) + n_sp*(log(n_sp) - log(T));

                ll2_model = ll_model./log(2);
                ll2_null  = ll_null./log(2);
                bitsPerSpike_test(ll,cc,hh,gg) = (ll2_model-ll2_null)./n_sp;
                
                %% pseudo-R^2
                
                pseudo_R2_test(ll,cc,hh,gg) = 1-(-n_sp - ll_model)./(-n_sp - ll_null);
%                 pseudo_R2_test(ll,cc,hh,gg) = 1-( ll_model)./(ll_null);
                
                %% simulate GLMs 
                stim = X_c*glm_k+glm_b;
                startSpikes = psth_hh(1:TT_start,ll,cc);
                
                [psth_c,ss] = kcSimGLMs(glm_h_filt,dt, 0, stim',startSpikes', NS, devNum);
                psth_glm(:,ll,cc,hh,gg) = psth_c;

                %% compare PSTH R^2
                
                psthSmooth_glm = conv2(psth_c(:),ff_smooth,'same');
                
                [R2_test(ll,cc,hh,gg),R_test(ll,cc,hh,gg)] = getR2(psthSmooth_hh(T_idx,ll,cc),psthSmooth_glm(T_idx));

                
                %% compare coherence
                    
                Fs = 1e3;

                [wcoh,wcs,pp] = wcoherence(psth_hh(:,ll,cc),psth_glm(:,ll,cc,hh,gg),Fs);
                mcohs = mean(wcoh(pp<=freqCutoff,T_idx),2);
                mcss  = mean(   (angle(wcs(pp<=freqCutoff,T_idx))),2);
                mcs2s = mean(abs(angle(wcs(pp<=freqCutoff,T_idx))),2);
                mcs3s = mean(abs(     (wcs(pp<=freqCutoff,T_idx))),2);


                c_test(:,ll,cc,hh,gg)  = mcohs;
                p_test(:,ll,cc,hh,gg)  = mcss;
                p2_test(:,ll,cc,hh,gg) = mcs2s;
                p3_test(:,ll,cc,hh,gg) = mcs3s;

            end
            
            
            %%
            for cc_c = 1:NC
                for ll_c = 1:NL
                    if(isnan(R2_test(ll_c,cc_c,1,1)))
                        continue;
                    end
                    subplot(3*NC,NL,ll_c+(cc_c-1)*NL);
                    semilogx(spkHistLengths,squeeze(R2_test(ll_c,cc_c,:,:)),'o-');
                    if(ll_c == 1)
                        xlabel('spk hist length (ms)');
                        ylabel('R^2');
                        if(cc_c == 1)
                            title('sine modulated')
                        elseif(cc_c == 2)
                            title('square modulated')
                        end
                    end
                        
                    subplot(3*NC,NL,ll_c+(cc_c-1)*NL+NC*NL);
                    semilogx(spkHistLengths,squeeze(pseudo_R2_test(ll_c,cc_c,:,:)),'o-');
                    if(ll_c == 1)
                        xlabel('spk hist length (ms)');
                        ylabel('pseudo R^2');
                        if(cc_c == 1)
                            title('sine modulated')
                        elseif(cc_c == 2)
                            title('square modulated')
                        end
                    end
                    
                    subplot(3*NC,NL,ll_c+(cc_c-1)*NL+2*NC*NL);
                    semilogx(spkHistLengths,squeeze(bitsPerSpike_test(ll_c,cc_c,:,:)),'o-');
                    if(ll_c == 1)
                        xlabel('spk hist length (ms)');
                        ylabel('bits per spike');
                        if(cc_c == 1)
                            title('sine modulated')
                        elseif(cc_c == 2)
                            title('square modulated')
                        end
                    end
                    
                    
                    
                end
            end
            drawnow;
        end
        
        
    end
end

if(~DEBUG)
    save(sprintf('Results/Lundstrom/glmFitMetrics_%d.mat',Period_const),'X2','T_idx','c_test','pseudo_R2_test','psth_glm','p_test','p2_test','p3_test','ll_test','R_test','R2_test','freqs','psth_hh','smooth_sig','bitsPerSpike_test','-v7.3');
end