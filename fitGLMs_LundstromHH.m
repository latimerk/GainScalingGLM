rng(08042018);

load("Results/Lundstrom/LundstromSims.mat");
addpath GLMtools;
addpath kGLM;
addGLMPaths


N_GPUs = 4;

dt = 1e-3;
N_P  = length(StimPeriods);
N_L  = length(StimLevels);


[stimBasis,stimBasis_0] = getStimBasis();

[~,sb] = getSpkHistBasis((16e3-1)/1e3,25);
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
    spkHistBases{H}   = orth(spkHistBases_0{H});
    
    N_Bspks(H)  = size(spkHistBases{H},2);
end

N_Bstim = size(stimBasis,2);
T_range = max(T_range_downsampled,2);

% XX = conv2(X1,stimBasis);
% XX = XX(1:length(X1),:);

TT = size(X1,1);
TT_fit = T_range(2) - T_range(1) + 1;




N_T = 6;
N_SPK = max(N_Bspks);
glm_k_stim = nan(N_Bstim,N_P+1,N_L,N_H, N_T);
glm_h_spk  = nan(N_SPK,N_P+1,N_L,N_H, N_T);
glm_b      = nan(N_P+1,N_L,N_H, N_T);
glm_ll     = nan(N_L,N_P,N_P+1,N_L,N_H, N_T);



simTypes = cell(N_T*2,1);

opts = optimoptions("fminunc","display","off","gradobj","on","hessian","on","maxiter",200,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu = optimoptions("fminunc","display","off","gradobj","on","hessian","on","maxiter",1000,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu2 = optimoptions("fminunc","display","off","gradobj","on","hessian","off","maxiter",200,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','quasi-newton');


sqGainFunc2   = @(t,Period,Height) round(sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
sineGainFunc2 = @(t,Period,Height)      (sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;


%%
T_start = 1;
L_start = 1;
H_start = 1;


%%
XX     = zeros(TT,N_P,N_Bstim,N_L);
XX_all = zeros(TT_fit,N_P,N_Bstim);
XX2    = zeros(TT_fit,N_P,N_Bstim,N_L);
for T = T_start:(N_T)
    %%
    fprintf('T = %d\n',T);
    switch T
        case 1
            sts_c = sts_sq_2Na;
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_2Na";
            continue;
            
        case 2
            sts_c = sts_sq_3AHP;
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_3AHP";
            
        case 3
            sts_c = sts_sq_Default;
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_Default";
            continue;
            
        case 4
            sts_c = sts_sine_2Na;
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_2Na";
            continue;
            
        case 5
            sts_c = sts_sine_3AHP;
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_3AHP";
            
        case 6
            sts_c = sts_sine_Default;
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_Default";
            continue;
    end
    
    

    for P_fit = 1:N_P
        for L_fit = 1:N_L
            GG = gainFunc(1:TT,StimPeriods(P_fit),StimLevels(L_fit));
            X_c = conv2(X1.*GG,stimBasis);
            XX(:,P_fit,:,L_fit)       = X_c(1:size(X1,1),:);
            
            XX2(:,P_fit,:,L_fit) = XX(T_range(1):T_range(2),P_fit,:,L_fit);
        end
    end
    clear GG
    clear X_c
    
    %%

    
    for L_fit = L_start:N_L
        fprintf('  L_fit = %d\n',L_fit);
        if(T <= N_T)
            XX_all(:,:,:) = XX2(:,:,:,L_fit);
        end
        %%
        for H = H_start:N_H
            fprintf('    H = %d\n',H);
            %%
            
            Y_all      = zeros(TT_fit,N_P);
            Y_hist_all = zeros(TT_fit,N_P,N_Bspks(H));
            spkHist_c = spkHistBases{H};
            N_Bspks_c = N_Bspks(H);
            %% fit GLMs to individual levels
            w_fits = zeros(N_Bstim+N_Bspks(H)+1,N_P+1);
            parfor P_fit = 1:N_P
                fprintf('      P_fit = %d\n',P_fit);
                %% setup spike hist
                T_range_all = (1:TT_fit) + TT_fit*(P_fit-1);

                Y = zeros(TT,1);
                Y(sts_c{P_fit,L_fit}) = 1;

                Y_hist = conv2(Y,spkHist_c);
                Y_hist_all(:,P_fit,:) = Y_hist((T_range(1):T_range(2))-1,:); %#ok<PFBNS>
                Y = Y(T_range(1):T_range(2));

                Y_all(:,P_fit) = Y;

                %% setup GLM
                X_glm = [squeeze(XX_all(:,P_fit,:)) squeeze(Y_hist_all(:,P_fit,:)) ones(TT_fit,1)];
                %Xy = X_glm'*Y;
                %normalizer = log(dt)*sum(Y)-sum(gammaln(Y+1));


                X_gpu = kcArrayToGPU((X_glm),mod(P_fit,N_GPUs));
                Y_gpu = kcArrayToGPU((Y),mod(P_fit,N_GPUs));
                C_gpu = kcArrayToGPU((zeros(size(X_glm,1),size(X_glm,2)+7)),mod(P_fit,N_GPUs));

                nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,0,X_gpu,Y_gpu,C_gpu);
%                 nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,1,X_gpu,Y_gpu,C_gpu);

                %% fit GLM
                %nllFunc = @(w) glmNll_PoissonExp(w,X_glm,Y,dt,Xy,normalizer);
                w_init = zeros(size(X_glm,2),1);


                w_fit_0 = fminunc(nllFuncGPU,w_init,opts_gpu); 
                kcFreeGPUArray(X_gpu);
                kcFreeGPUArray(Y_gpu);
                kcFreeGPUArray(C_gpu);


                w_fits(:,P_fit) = w_fit_0;%fminunc(nllFunc,w_fit_0,opts); %in the unlikely case single precision is a problem, finishes fit on CPU with double precision


            end
            for P_fit = 1:N_P
                glm_k_stim(:,P_fit,L_fit,H,T) = w_fits(1:N_Bstim,P_fit);
                hh = [w_fits((1:N_Bspks_c)+N_Bstim,P_fit);nan(N_SPK-N_Bspks_c,1)];
                glm_h_spk(:,P_fit,L_fit,H,T) = hh;
                glm_b(P_fit,L_fit,H,T) = w_fits(1+N_Bspks_c+N_Bstim,P_fit);
            end

            %% fit GLM to all periods
            X_glm = [reshape(XX_all,[],N_Bstim) reshape(Y_hist_all,[],N_Bspks_c) ones(numel(Y_all),1)];

            Y_all_flat = Y_all(:);
            %Xy = X_glm'*Y_all_flat;
            %normalizer = log(dt)*sum(Y_all(:))-sum(gammaln(Y_all(:)+1));
            w_init = zeros(size(X_glm,2),1);

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

            nllFuncGPU = @(w)kcGlmNLL_Multi(w,dt,0,X_gpu,Y_gpu,C_gpu);


            %nllFunc = @(w) glmNll_PoissonExp(w,X_glm,Y_all_flat,dt,Xy,normalizer);  
            w_fit_0 = fminunc(nllFuncGPU,w_init,opts_gpu);  
            w_fits(:,N_P+1) = w_fit_0;%fminunc(nllFunc,w_fit_0,opts); %in the unlikely case single precision is a problem, finishes fit on CPU with double precision

            P_fit = N_P + 1;
            hh = [w_fits((1:N_Bspks_c)+N_Bstim,P_fit);nan(N_SPK-N_Bspks_c,1)];

            glm_k_stim(:,P_fit,L_fit,H,T) = w_fits(1:N_Bstim,N_P+1);
            glm_h_spk (:,P_fit,L_fit,H,T) = hh;
            glm_b(P_fit,L_fit,H,T) = w_fits(1+N_Bspks(H)+N_Bstim,N_P+1);


            for gg = 1:N_GPUs
                kcFreeGPUArray(X_gpu{gg});
                kcFreeGPUArray(Y_gpu{gg});
                kcFreeGPUArray(C_gpu{gg});
            end

            clear X_glm;
            clear Y_all_flat;

            %% get LL of fits for each period & level & simualate

            X_glm = ones(T_range(2)-T_range(1)+1,size(XX,3)+size(Y_hist_all,3)+1);
            for P_test = 1:N_P
                fprintf('     P_test = %d\n',P_test);
                X_glm(:,size(XX,3)+(1:N_Bspks_c)) = squeeze(Y_hist_all(:,P_test,:));
                for L_test = 1:N_L
                    fprintf('        L_test = %d\n',L_test);
                    %%
                
                    X_glm(:,1:size(XX,3)) = squeeze(XX2(:,P_test,:,L_test)) ;
                    Y     = (Y_all(:,P_test));

                    Xy = X_glm'*Y;
                    normalizer = log(dt)*sum(Y)-sum(gammaln(Y+1));

                    
                    
                    eX = exp(min(700,(X_glm*w_fits) + log(dt)));

                    f = -(-Xy'*w_fits + sum(eX) - normalizer);
                    glm_ll(L_test,P_test,:,L_fit,H,T) = f;

                end
            end
            clear X_glm;
            
            %%
            clear Y_hist_all;
            if(mod(H,4) == 0 || H == N_H)
                save("Results/Lundstrom/LundstromGLMs_meta.mat","-v7.3","T","H","L_fit","simTypes","glm_k_stim","glm_h_spk","glm_b","glm_ll","stimBasis","spkHistBases","spkHistLengths","spkHistBases_0");
            end
        end
        H_start = 1;
    end
    L_start = 1;

end

