rng(08042018);

load("Results/Lundstrom/LundstromSims.mat");
addpath GLMtools;
addpath ~/gitCode/GLM/
addGLMPaths


load('Results/Lundstrom/LundstromGLMs_meta.mat');
N_GPUs = 4;

dt = 1e-3;
N_P  = length(StimPeriods);
N_L  = length(StimLevels);
N_H  = length(spkHistLengths);

N_Bspks = zeros(N_H,1);
for H = 1:N_H
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



simTypes = cell(N_T*2,1);

opts = optimoptions("fminunc","display","final","gradobj","on","hessian","on","maxiter",200,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu = optimoptions("fminunc","display","final","gradobj","on","hessian","on","maxiter",1000,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','trust-region');
opts_gpu2 = optimoptions("fminunc","display","off","gradobj","on","hessian","off","maxiter",200,"StepTolerance",1e-6,"FunctionTolerance",1e-6,'algorithm','quasi-newton');


sqGainFunc2   = @(t,Period,Height) round(sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
sineGainFunc2 = @(t,Period,Height)      (sin( (((t(:)-BUFFER/downsampleRate)*h*downsampleRate)/1000.0 *2*pi)./Period).*0.5+0.5).*(Height(:)'-1) + 1;
%%



allT = 1:(N_T*2);

%%
XX     = zeros(TT,N_P,N_Bstim,N_L);
% TS_to_fit = 1:12;
for T = TS_to_fit
    fprintf('T = %d\n',T);
    if(exist(sprintf('Results/Lundstrom/LundstromGLMs_part%d.mat',T),'file'))
        fprintf('Loading partial results...\n');
        load(sprintf('Results/Lundstrom/LundstromGLMs_part%d.mat',T));
        
        loadedSts = true;
        L_start = L_fit;
        H_start = H  +  1;
    else
        
        L_start = 1;
        H_start = 1;
        loadedSts = false;
    end
    
    %%
    switch T
        case 1
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_2Na";
            
            if(~loadedSts)
                sts_glm_sq_2Na     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 2
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_3AHP";
            
            if(~loadedSts)
                sts_glm_sq_3AHP     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
        case 3
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_Default";
            
            if(~loadedSts)
                sts_glm_sq_Default     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 4
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_2Na";
            
            if(~loadedSts)
                sts_glm_sine_2Na     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 5
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_3AHP";
            
            if(~loadedSts)
                sts_glm_sine_3AHP     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
        case 6
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_Default";
            
            
            if(~loadedSts)
                sts_glm_sine_Default     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 7
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_sqTrain_2Na";
            
            if(~loadedSts)
                sts_glm_sine_sqTrain_2Na     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 8
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_sqTrain_3AHP";
            
            if(~loadedSts)
                sts_glm_sine_sqTrain_3AHP     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
        case 9
            gainFunc = sineGainFunc2;
            simTypes{T} = "sine_sqTrain_Default";
            
            if(~loadedSts)
                sts_glm_sine_sqTrain_Default     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 10
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_sineTrain_2Na";
            
            if(~loadedSts)
                sts_glm_sq_sineTrain_2Na     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
        case 11
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_sineTrain_3AHP";
            
            if(~loadedSts)
                sts_glm_sq_sineTrain_3AHP     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
        case 12
            gainFunc = sqGainFunc2;
            simTypes{T} = "sq_sineTrain_Default";
            if(~loadedSts)
                sts_glm_sq_sineTrain_Default     = cell(N_L,N_P,N_P+1,N_L,N_H);
            end
            continue;
            
    end
    
    

    for P_fit = 1:N_P
        for L_fit = 1:N_L
            GG = gainFunc(1:TT,StimPeriods(P_fit),StimLevels(L_fit));
            X_c = conv2(X1.*GG,stimBasis);
            XX(:,P_fit,:,L_fit)       = X_c(1:size(X1,1),:);
        end
    end
    clear GG
    clear X_c
    
    %%
    
    for L_fit = L_start:N_L
        fprintf('  L_fit = %d\n',L_fit);
        %%
        for H = H_start:N_H
            fprintf('    HH = %d\n',H);
            %%
            
            
            spkHist_c = spkHistBases{H};
            N_Bspks_c = N_Bspks(H);

            %% get LL of fits for each period & level & simualate
            %XX4 = XX(1:10e3,:,:,:);
            tic;
            XX_c = zeros(size(XX,1),N_L,N_P+1);
            XX2_c = zeros(size(XX,1),size(stimBasis,2),N_L);
            for P_test = 1:N_P 
                fprintf('      Starting sim round %d.\n',P_test);
                
                sts_sim = cell(N_L,1,N_P+1);
                
                gk = glm_k_stim(:,:,L_fit,H,mod(T-1,N_T)+1);
                
                XX2_c = squeeze(XX(:,P_test,:,:));
                for L_test = 1:N_L
                    XX_c(:,L_test,:) = XX2_c(:,:,L_test)*gk;
                end
                %%
                
                parfor P_fit = 1:(N_P+1)

                    %% run simulation
                    
                    h_filt      = spkHist_c*glm_h_spk( 1:N_Bspks_c,P_fit,L_fit,H,mod(T-1,N_T)+1); 
                    b     = glm_b(P_fit,L_fit,H,mod(T-1,N_T)+1);
                    
                    %bb = tic;
                    sts_c = simGLMs_cpu(h_filt,dt, 0, 1, XX_c(:,:,P_fit)'+b);
                    
                    sts_sim(:,1,P_fit) = sts_c;
                    %s = toc(bb);
                    %fprintf('seconds per GLM = %.2f\n',s/size(XX_c,2));
                    
                end
                
                %%
                fprintf('      Finishing sim round %d.\n',P_test);
                switch T
                    case 1
                        sts_glm_sq_2Na(:,P_test,:,L_fit,H) = sts_sim;
                    case 2
                        sts_glm_sq_3AHP(:,P_test,:,L_fit,H) = sts_sim;
                    case 3
                        sts_glm_sq_Default(:,P_test,:,L_fit,H) = sts_sim;
                    case 4
                        sts_glm_sine_2Na(:,P_test,:,L_fit,H) = sts_sim;
                    case 5
                        sts_glm_sine_3AHP(:,P_test,:,L_fit,H) = sts_sim;
                    case 6
                        sts_glm_sine_Default(:,P_test,:,L_fit,H) = sts_sim;


                    case 7
                        sts_glm_sine_sqTrain_2Na(:,P_test,:,L_fit,H) = sts_sim;
                    case 8
                        sts_glm_sine_sqTrain_3AHP(:,P_test,:,L_fit,H) = sts_sim;
                    case 9
                        sts_glm_sine_sqTrain_Default(:,P_test,:,L_fit,H) = sts_sim;
                    case 10
                        sts_glm_sq_sineTrain_2Na(:,P_test,:,L_fit,H) = sts_sim;
                    case 11
                        sts_glm_sq_sineTrain_3AHP(:,P_test,:,L_fit,H) = sts_sim;
                    case 12
                        sts_glm_sq_sineTrain_Default(:,P_test,:,L_fit,H) = sts_sim;
                        
                end


                clear sts_sim
            end
            
            toc
            
            %%
            clear Y_hist_all;
            if(mod(H,4) == 0 || H == N_H)
                
                switch T
                    case 1
                        save("Results/Lundstrom/LundstromGLMs_part1.mat","-v7.3","sts_glm_sq_2Na","T","L_fit","H");
                    case 2
                        save("Results/Lundstrom/LundstromGLMs_part2.mat","-v7.3","sts_glm_sq_3AHP","T","L_fit","H");
                    case 3
                        save("Results/Lundstrom/LundstromGLMs_part3.mat","-v7.3","sts_glm_sq_Default","T","L_fit","H");
                    case 4
                        save("Results/Lundstrom/LundstromGLMs_part4.mat","-v7.3","sts_glm_sine_2Na","T","L_fit","H");
                    case 5
                        save("Results/Lundstrom/LundstromGLMs_part5.mat","-v7.3","sts_glm_sine_3AHP","T","L_fit","H");
                    case 6
                        save("Results/Lundstrom/LundstromGLMs_part6.mat","-v7.3","sts_glm_sine_Default","T","L_fit","H");
                    case 7
                        save("Results/Lundstrom/LundstromGLMs_part7.mat","-v7.3","sts_glm_sine_sqTrain_2Na","T","L_fit","H");
                    case 8
                        save("Results/Lundstrom/LundstromGLMs_part8.mat","-v7.3","sts_glm_sine_sqTrain_3AHP","T","L_fit","H");
                    case 9
                        save("Results/Lundstrom/LundstromGLMs_part9.mat","-v7.3","sts_glm_sine_sqTrain_Default","T","L_fit","H");
                    case 10
                        save("Results/Lundstrom/LundstromGLMs_part10.mat","-v7.3","sts_glm_sq_sineTrain_2Na","T","L_fit","H");
                    case 11
                        save("Results/Lundstrom/LundstromGLMs_part11.mat","-v7.3","sts_glm_sq_sineTrain_3AHP","T","L_fit","H");
                    case 12
                        save("Results/Lundstrom/LundstromGLMs_part12.mat","-v7.3","sts_glm_sq_sineTrain_Default","T","L_fit","H");
                end
            end
        end
        H_start = 1;
    end
    L_start = 1;
    
    switch T
        case 1
            clear sts_glm_sq_2Na
        case 2
            clear sts_glm_sq_3AHP
        case 3
            clear sts_glm_sq_Default
        case 4
            clear sts_glm_sine_2Na
        case 5
            clear sts_glm_sine_3AHP
        case 6
            clear sts_glm_sine_Default
        case 7
            clear sts_glm_sine_sqTrain_2Na
        case 8
            clear sts_glm_sine_sqTrain_3AHP
        case 9
            clear sts_glm_sine_sqTrain_Default;
        case 10
            clear sts_glm_sq_sineTrain_2Na;
        case 11
            clear sts_glm_sq_sineTrain_3AHP;
        case 12
            clear sts_glm_sq_sineTrain_Default;
    end
end

