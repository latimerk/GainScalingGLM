alphas_all = [0 1:0.5:5 10];

alphas = alphas_all(SIM_ALPHAS);
AA = length(alphas);
addpath Utils
addpath EMD/

for alphaCtr = 1:AA
    alpha = alphas(alphaCtr);
    if(alpha <= 0)
        saveDir = sprintf('Results/Mease/');
    else
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    end

    load(sprintf('%s/Mease_gainInfo.mat',saveDir));
    load(sprintf('%s/MeaseGLMs.mat',saveDir));

    load('Results/Mease/MeaseSims.mat', 'stim_mu','xs','frFunctions', 'G_mat_info_paper')

    spkHistLengths(1) = 10;

    fr_0 = squeeze(frFunctions(xs==0,:,:));
    fprintf('calculating distance metrics...\n');
    H = length(spkHistLengths);

    gainDistance = struct();

    gainDistance.hh.JS = nan(15,15,3);
    gainDistance.hh.EM = nan(15,15,3);

    gainDistance.glm.JS = nan(15,15,3,H,5);
    gainDistance.glm.EM = nan(15,15,3,H,5);
    
    gainDistance.alpha = alpha;

    for nn = 1:15
        fprintf('  nn = %d / 15\n',nn);
        for kk = 1:15
            fprintf('    kk = %d / 15\n',kk);

            hh_base = fits_gain.hh.p_spk(:,nn,kk,1);
            if(fr_0(nn,kk) > 1)
                continue;
            end

            for cc = 1:3
                hh_comp = fits_gain.hh.p_spk(:,nn,kk,cc+1);
                dx = fits_gain.bins(2)-fits_gain.bins(1);
                gainDistance.hh.JS(nn,kk,cc) = getJS(hh_base,hh_comp);
                gainDistance.hh.EM(nn,kk,cc) = getWD(hh_base,hh_comp,fits_gain.bins);


            end

            
            gd_JS = nan(3,H,5);
            gd_EM = nan(3,H,5);
            parfor hh = 1:H
                for tt = 1:5
                    glm_base = fits_gain.glm.p_spk(:,nn,kk,1,tt,hh);
                    for cc = 1:3
                        glm_comp = fits_gain.glm.p_spk(:,nn,kk,cc+1,tt,hh);


%                         gainDistance.glm.JS(nn,kk,cc,hh,tt) = getJS(glm_base,glm_comp);
%                         gainDistance.glm.EM(nn,kk,cc,hh,tt) = getWD(glm_base,glm_comp,fits_gain.bins);

                        dx = fits_gain.bins(2)-fits_gain.bins(1);
                        gd_JS(cc,hh,tt) = getJS(glm_base,glm_comp);
                        gd_EM(cc,hh,tt) = getWD(glm_base,glm_comp,fits_gain.bins);
                    end

                end
            end
            
            gainDistance.glm.JS(nn,kk,:,:,:) = gd_JS;
            gainDistance.glm.EM(nn,kk,:,:,:) = gd_EM;
        end
    end
    fprintf('  saving...\n');
    distanceFile = sprintf('%s/gainDistance.mat',saveDir);
    save(distanceFile,'gainDistance');
    fprintf('done.\n');
end