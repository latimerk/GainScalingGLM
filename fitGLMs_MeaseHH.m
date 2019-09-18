

load("Results/Mease/MeaseSims.mat");
addpath GLMtools;
addpath kGLM;
addGLMPaths;

[~,ii] = min(abs(xs-0));
fr0 = squeeze(frFunctions(ii,:,:));

alphas_all = [0 1:0.5:5 10];

alphas = alphas_all(SIM_ALPHAS);
AA = length(alphas);

if(~exist('IDX','var'))
    IDX = randi(100);
end


N_GPUs = 4;

T_range_downsampled = max(T_range_downsampled,2);



dt = 1e-3;
N_Na = size(sts,1);
N_K  = size(sts,2);
N_L  = length(StimLevels);

stimBasis = getStimBasis();
[~,sb] = getSpkHistBasis((150-1)/1e3,15);
[~,spkHistLengths] = max(sb);
spkHistLengths = spkHistLengths(5:end);
spkHistLengths(1) = 10;

N_H  = length(spkHistLengths);
spkHistBases   = cell(N_H,1);
spkHistBases_0   = cell(N_H,1);
N_Bspks = zeros(N_H,1);


for H = 1:N_H
    %[spkHistBases{H},spkHistBases_0{H}] = getSpkHistBasis(spkHistLengths(H)*1e-3);
    
    hh = sb(:,1:H+4);
    
    th = sum(abs(hh),2);
    tt = find(th>0,1,'last');
    if(isempty(tt))
        tt = size(sb,1);
    end
    
    spkHistBases_0{H} = hh(1:tt,:);
    spkHistBases{H}   = [orth(spkHistBases_0{H}(:,1:5)) orth(spkHistBases_0{H}(:,6:end))];
    
    N_Bspks(H)  = size(spkHistBases{H},2);
end

N_Bstim = size(stimBasis,2);

XX = conv2(X1,stimBasis);
XX = XX(1:length(X1),:);

TT = size(XX,1);
TT_fit = T_range_downsampled(2) - T_range_downsampled(1) + 1;

XX_all = zeros(TT_fit*N_L,N_Bstim);
XX_cs = zeros(TT_fit,N_Bstim,N_L);
for L = 1:N_L
    XX_cs(:,:,L) = XX(T_range_downsampled(1):T_range_downsampled(2),:)*StimLevels(L);
    XX_all((1:TT_fit)+(L-1)*TT_fit,:) = XX_cs(:,:,L);
end

glm_k_stims = nan(N_Bstim,N_Na,N_K,N_L+1,N_H,AA);
glm_h_spks  = nan(max(N_Bspks),N_Na,N_K,N_L+1,N_H,AA);
glm_bs      = nan(N_Na,N_K,N_L+1,N_H,AA);
glm_lls     = nan(N_L,N_Na,N_K,N_L+1,N_H,AA);

opts = optimoptions(     "fminunc","display","off","gradobj","on","hessian","on" ,"maxiter",200 ,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu = optimoptions( "fminunc","display","off","gradobj","on","hessian","on" ,"maxiter",1000,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu2 = optimoptions("fminunc","display","off","gradobj","on","hessian","off","maxiter",200 ,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','quasi-newton');

G_Na_init = 1;
G_K_init = 1;
H_init = 1;
G_K = 0;

%%
if(exist(sprintf("Results/Mease/MeaseGLMs_multi_part_%d.mat",IDX),'file'))
    load(sprintf("Results/Mease/MeaseGLMs_multi_part_%d.mat",IDX));
    G_Na_init = G_Na;
    G_K_init = G_K+1;
end
%%



for G_Na = G_Na_init:N_Na
    fprintf('Na = %d / %d\n',G_Na,N_Na);
    for G_K = G_K_init:N_K
        fprintf('   K = %d / %d\n',G_K,N_K);

        tic; 
        if(isnan(stim_mu(G_Na,G_K)) || length(sts{G_Na,G_K,1}) < 10 || fr0(G_Na,G_K) > 1)
            continue;
        end
        
        %%
        
        Y_0_all      = zeros(TT,N_L,1);
        Y_all      = zeros(TT_fit,N_L,1);
        for L = 1:N_L
            Y_0_all(sts{G_Na,G_K,L},L) = 1;
            Y_all(:,L) = Y_0_all(T_range_downsampled(1):T_range_downsampled(2),L);
        end
                
        for H = H_init:N_H
            
            %%
            fprintf('      H = %d / %d\n',H,N_H);
            spkHistBasis_c = spkHistBases{H};
            N_Bspks_c = N_Bspks(H);
            Y_hist_all = zeros(TT_fit*N_L,N_Bspks_c);
            Y_hists    = zeros(TT_fit    ,N_Bspks_c,N_L);
            w_fits = zeros(N_Bstim+N_Bspks_c+1,N_L+1,AA);

            
            for L = 1:N_L
                Y_hist = conv2(Y_0_all(:,L),spkHistBasis_c);
                Y_hist = Y_hist((T_range_downsampled(1):T_range_downsampled(2))-1,:);
                Y_hist_all((1:TT_fit)+(L-1)*TT_fit,:) = Y_hist;
                Y_hists(:,:,L) = Y_hist;
            end
            N_GPUs_c = 1;
            
            %%
            for L = 1:N_L
                fprintf('        L = %d / (%d + 1)\n',L,N_L);
                alphas_c = alphas;
                %% setup GLM

                
                Y = Y_all(:,L);

                X_glm = [XX_cs(:,:,L) Y_hists(:,:,L) ones(TT_fit,1)];
                Xy = X_glm'*Y;

                normalizer = log(dt)*sum(Y)-sum(gammaln(Y+1));

                %% fit GLM

                devNums = mod((0:(N_GPUs_c-1)) + (L-1)*N_GPUs_c,N_GPUs);

                NperGPU = round(size(X_glm,1)/N_GPUs_c);
                X_gpu = cell(N_GPUs_c,1);
                Y_gpu = cell(N_GPUs_c,1);
                C_gpu = cell(N_GPUs_c,1);
                for gg = 1:N_GPUs_c
                    kcResetDevice(devNums(gg));
                    endT    = gg*NperGPU;
                    if(gg == N_GPUs_c)
                        endT = size(X_glm,1);
                    end
                    tts = ((gg-1)*NperGPU+1):endT;
                    X_gpu{gg} = kcArrayToGPU(X_glm(tts,:),devNums(gg));
                    Y_gpu{gg} = kcArrayToGPU(Y(tts),devNums(gg));
                    C_gpu{gg} = kcArrayToGPU(zeros(length(tts),size(X_glm,2)+7),devNums(gg));
                end
                %%
                w_0s = zeros(size(X_glm,2),AA);
                
                for alphaCtr = AA:-1:1
                    %%

                    alpha = alphas_c(alphaCtr);


                    if(alpha > 0)
                        if(alpha == 1)
                            GLMtype = 4;%Poisson: 4, Ber: 2
                            nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,GLMtype,X_gpu,Y_gpu,C_gpu);
                            %nllFunc = @(w) glmNll_PoissonSoftRec(w,X_glm,Y,dt,normalizer);
                        else
                            GLMtype = 5;%Poisson:5, Ber: 3
                            nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,GLMtype,X_gpu,Y_gpu,C_gpu,alpha);
                            %nllFunc = @(w) glmNll_PoissonSoftPow(w,X_glm,Y,alpha,dt,normalizer);
                        end
                    else
                        GLMtype = 0;%Poisson: 0, Ber: 1
                        nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,GLMtype,X_gpu,Y_gpu,C_gpu);
                    end

                    if(alpha <= 0 || alphaCtr == AA)
                        w_init = zeros(size(X_glm,2),1);
                    else
                        w_init = w_0s(:,alphaCtr+1);
                    end


                    w_0s(:,alphaCtr) = fminunc(nllFuncGPU,w_init,opts_gpu);
                    
                end
                
                w_fits(:,L,:) = w_0s;
                %%
                for gg = 1:N_GPUs_c
                    kcFreeGPUArray(X_gpu{gg});
                    kcFreeGPUArray(Y_gpu{gg});
                    kcFreeGPUArray(C_gpu{gg});
                    kcResetDevice(devNums(gg));
                end

            end
            for alphaCtr = 1:AA
                for L = 1:N_L
                    glm_k_stims(:         ,G_Na,G_K,L,H,alphaCtr) = w_fits(1:N_Bstim            ,L,alphaCtr);
                    glm_h_spks(1:N_Bspks_c,G_Na,G_K,L,H,alphaCtr) = w_fits((1:N_Bspks_c)+N_Bstim,L,alphaCtr);
                    glm_bs(                G_Na,G_K,L,H,alphaCtr) = w_fits( 1+N_Bspks_c +N_Bstim,L,alphaCtr);

                end
            end

            %% fit GLM to all levels
            fprintf('        L = (%d + 1) / (%d + 1)\n',N_L,N_L);
            X_glm = [XX_all Y_hist_all ones(TT_fit*N_L,1)];
            Y_all_flat = Y_all(:);
            %Xy = X_glm'*Y_all_flat;
            %normalizer = log(dt)*sum(Y_all(:))-sum(gammaln(Y_all(:)+1));


            X_gpu = cell(N_GPUs,1);
            Y_gpu = cell(N_GPUs,1);
            C_gpu = cell(N_GPUs,1);
            T_per = ceil(size(X_glm,1)/N_GPUs);
            for gg = 1:N_GPUs
                tts = (1:T_per) + (gg-1)*T_per;
                if(gg == N_GPUs)
                    tts = tts(1):size(X_glm,1);
                end
                X_gpu{gg} = kcArrayToGPU((X_glm(tts,:)),gg-1);
                Y_gpu{gg} = kcArrayToGPU((Y_all_flat(tts)),gg-1);
                C_gpu{gg} = kcArrayToGPU((zeros(length(tts),size(X_glm,2)+7)),gg-1);
            end

            for alphaCtr = 1:AA
                %%
                alpha = alphas(alphaCtr);
                if(alpha > 0)
                    if(alpha == 1)
                        nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,4,X_gpu,Y_gpu,C_gpu);
                    else
                        nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,5,X_gpu,Y_gpu,C_gpu,alpha);
                    end
                else
                    nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,0,X_gpu,Y_gpu,C_gpu);
                end

                w_init = zeros(size(X_glm,2),1);
                w_fits(:,N_L+1,alphaCtr) = fminunc(nllFuncGPU,w_init,opts_gpu); 

                glm_k_stims(:,         G_Na,G_K,N_L+1,H,alphaCtr) = w_fits(1:N_Bstim            ,N_L+1,alphaCtr);
                glm_h_spks(1:N_Bspks_c,G_Na,G_K,N_L+1,H,alphaCtr) = w_fits((1:N_Bspks_c)+N_Bstim,N_L+1,alphaCtr);
                glm_bs(                G_Na,G_K,N_L+1,H,alphaCtr) = w_fits(1+N_Bspks_c+N_Bstim  ,N_L+1,alphaCtr);
            end
            
            for gg = 1:N_GPUs
                kcFreeGPUArray(X_gpu{gg});
                kcFreeGPUArray(Y_gpu{gg});
                kcFreeGPUArray(C_gpu{gg});
                kcResetDevice(gg-1);
            end


            clear Y_all_flat;
        end

        %%
        H_init = 1;
        tm = toc;

        fprintf('\n\ntime to fit one G_K/G_Na set: %.2f\n\n\n',tm);

        save(sprintf("Results/Mease/MeaseGLMs_multi_part_%d.mat",IDX),"-v7.3","G_Na","G_K","glm_k_stims","glm_h_spks","glm_bs","glm_lls","stimBasis","spkHistBases","spkHistLengths","StimLevels");
        pause(30);
    end
    pause(60);
    G_K_init = 1;
end


%%
fprintf('saving final results...\n');
for alphaCtr = 1:AA
    alpha = alphas(alphaCtr);
    if(alpha > 0)
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    else
        saveDir = sprintf('Results/Mease/');
    end
    if(~isfolder(saveDir))
        mkdir(saveDir);
    end
    

    glm_k_stim = glm_k_stims(:,:,:,:,:,alphaCtr);
    glm_h_spk  = glm_h_spks(:,:,:,:,:,alphaCtr);
    glm_b      = glm_bs(:,:,:,:,alphaCtr);
    glm_ll     = glm_lls(:,:,:,:,:,alphaCtr);
    save(sprintf("%sMeaseGLMs.mat",saveDir),"-v7.3","alpha","glm_k_stim","glm_h_spk","glm_b","glm_ll","stimBasis","spkHistBases","spkHistLengths","StimLevels");
end
fprintf('done\n');