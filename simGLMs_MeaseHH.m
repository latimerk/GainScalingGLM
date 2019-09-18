rng(08032018);
addpath GLMtools;
poolobj = gcp('nocreate'); delete(poolobj);
parpool(5);

NO_OVERRIDE = false;

alphas_all = [0 1:0.5:5 10];

alphas = alphas_all(SIM_ALPHAS);

for alphaCtr = 1:length(alphas)

    alpha = alphas(alphaCtr);

    if(alpha > 0)
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    else
        saveDir = sprintf('Results/Mease/');
    end


    load(sprintf("%sMeaseGLMs.mat",saveDir));
    load("Results/Mease/MeaseSims.mat","X1","T_range_downsampled","xs","frFunctions");


    N_Na = size(glm_k_stim,2);
    N_K  = size(glm_k_stim,3);
    N_L  = length(StimLevels);
    N_H  = length(spkHistLengths);

    dt = 1e-3;

    N_Bspks = zeros(N_H,1);
    for H = 1:N_H
        N_Bspks(H)  = size(spkHistBases{H},2);
    end



    XX = conv2(X1,stimBasis);
    XX = XX(1:length(X1),:);
    TT = size(XX,1);
    TT_fit = T_range_downsampled(2) - T_range_downsampled(1) + 1;

    for G_Na = 1:N_Na
        fprintf('Na = %d / %d\n',G_Na,N_Na);
        for G_K = 1:N_K
            fprintf('   K = %d / %d\n',G_K,N_K);

            tic; 
            if(isnan(glm_k_stim(1,G_Na,G_K,1,1)))
                continue;
            end

            fName = sprintf("%sMeaseGLMs_spks_Na_%d_K_%d.mat",saveDir,G_Na,G_K);

            if(exist(fName,'file') && NO_OVERRIDE)
                fprintf('       results file found.\n');
                continue;
            end

            sts_glm = cell(N_L,N_L+1,N_H);                    
                    
            
            for H = 1:N_H
                fprintf('      H(sim) = %d / %d\n',H,N_H);
                sts_c = cell(N_L,N_L+1);
                spkHistBasis_c = spkHistBases{H};
                N_Bspks_c = N_Bspks(H);
                StimLevels_c = StimLevels(:);
                
                %%
                XXs = (XX*squeeze(glm_k_stim(:,G_Na,G_K,:,H)));
                
                parfor L_train = 1:(N_L+1)
                    %% run simulation

                    h_spk  = spkHistBasis_c*squeeze(glm_h_spk(1:N_Bspks_c,G_Na,G_K,L_train,H));

                    XX_c = StimLevels_c*(XXs(:,L_train))' + glm_b(G_Na,G_K,L_train,H);
                    
                    sts_c(:,L_train) = simGLMs_cpu(h_spk,dt, alpha, 1, XX_c);

                end

                sts_glm(:,:,H) = sts_c;
            end
                    
                    
                    
            tm = toc;

            fprintf('\n\ntime to sim one G_K/G_Na set: %.2f\n\n\n',tm);
            save(fName,"-v7.3","G_Na","G_K","sts_glm");
        end
    end
    fprintf('done.\n');
end