

%alpha = 1.5;

if(~exist('sts_0','var'))
    load('Results/Mease/MeaseSims.mat');
end
addpath Utils

alphas_all = [0 1:0.5:5 10];

alphas = alphas_all(SIM_ALPHAS);
AA = length(alphas);

for alphaCtr = 1:AA
    %%
    alpha = alphas(alphaCtr);

    if(alpha <= 0)
        saveDir = sprintf('Results/Mease/');
    else
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    end

    load(sprintf("%sMeaseGLMs.mat",saveDir))

    T_STA = 150;
    DT_STA = 0;


    N_Na = size(sts,1);
    N_K  = size(sts,2);
    N_L  = length(StimLevels);
    N_H  = length(spkHistLengths);

    % sta_glm = nan(T_STA,N_Na,N_K,N_L,  N_L+1,N_H);
    sta_hh = nan(T_STA,N_Na,N_K,N_L);
    sta_n_hh = nan(T_STA,N_Na,N_K,N_L);
    stc_hh = nan(T_STA,T_STA,N_Na,N_K,N_L);


    sta_glm = nan(T_STA,N_Na,N_K,N_L,N_L+1,N_H);
    sta_n_glm = nan(T_STA,N_Na,N_K,N_L,N_L+1,N_H);


    sta_avg_hh = nan(T_STA,N_Na,N_K,N_L);


    bins = (-5:0.1:5);

    s_spk_hh = nan(length(bins),N_Na,N_K,N_L);
    s_all_hh = nan(length(bins),N_Na,N_K,N_L);
    s_spk_avg_hh = nan(length(bins),N_Na,N_K,N_L,N_L);
    r_bar_hh = nan(N_Na,N_K,N_L);
    r_bar_glm = nan(N_Na,N_K,N_L,N_L+1,N_H);

    s_spk_glm = nan(length(bins),N_Na,N_K,N_L,N_L+1,N_H);
    s_all_glm = nan(length(bins),N_Na,N_K,N_L,N_L+1,N_H);
    k_spk_glm = nan(length(bins),N_Na,N_K,N_L,N_L+1,N_H);


    f_exp = exp(-(0:20)'./1);  f_exp = f_exp(:)./sum(f_exp);
    s_base = histc(X1,bins);


    isi_bins = 0:1:500;
    glmStats_isis = zeros(length(isi_bins),N_L,N_Na,N_K,N_L+1,N_H);
    hh_isis       = zeros(length(isi_bins),N_L,N_Na,N_K);

    ac_spk_hh = nan(N_Na,N_K,N_L);

    gain_bins1  = bins;
    gain_bins2  = 0:0.1:5;
%     if(alpha <= 0)
%         kh_spk_glm = zeros(length(gain_bins1),length(gain_bins2),N_Na,N_K,N_L,N_L+1,N_H,'single')-1;
%         kh_all_glm = zeros(length(gain_bins1),length(gain_bins2),N_Na,N_K,N_L,N_L+1,N_H,'single')-1;
%     else
        kh_spk_glm = [];
        kh_all_glm = [];
%     end


    TT = T_range_downsampled(2)-T_range_downsampled(1)+1;
    
    %%
    for nn = 1:N_Na
        fprintf('Na: %d / %d\n',nn,N_Na);
        for kk = 1:N_K
            X1_0 = X1;
            fprintf('  K: %d / %d\n',kk,N_K);
            %%
            %X_spk = cell(N_L,1);
            fName = sprintf("%sMeaseGLMs_spks_Na_%d_K_%d.mat",saveDir,nn,kk);
            if(~exist(fName,'file') )
                continue;
            else
                if(isnan(glm_b(nn,kk,1,1)))
                    fprintf('deleting invalid sim file for Na = %d, K = %d\n',nn,kk);
                    delete(fName);
                    continue;
                end
            end
            ff = load(fName,'sts_glm');
            sts_glm = ff.sts_glm;
            
            %%

            for N_L_test = 1:N_L
                fprintf('    N_L_test: %d / %d\n',N_L_test,N_L);
                spkHistBases_c = spkHistBases;
                sts_c = sts{nn,kk,N_L_test}(:);
                sts_c = sts_c(sts_c >= T_range_downsampled(1) & sts_c <= T_range_downsampled(2));

                dc = diff(sts_c);
                if(~isempty(dc))
                    hh_isis(:,N_L_test,nn,kk) =  histc(dc,isi_bins);
                else
                    hh_isis(:,N_L_test,nn,kk) = 0;
                end


                if(~isempty(sts_c))

                    N_spk = length(sts_c);
                    r_bar_hh(nn,kk,N_L_test) = N_spk./(T_range_downsampled(2)-T_range_downsampled(1) + 1)*1e3;

                    X_spk = X1_0(sts_c + DT_STA + (0:-1:(-T_STA+1)));
                    ss = mean(X_spk,1)';
                    sta_hh(:,nn,kk,N_L_test) = ss;
                    stc_hh(:,:,nn,kk,N_L_test) = cov(X_spk);

                    ss_n = ss./sqrt(sum(ss.^2));


                    sta_n_hh(:,nn,kk,N_L_test) = ss_n;

                    d = X_spk*ss_n;


                    s_spk_hh(:,nn,kk,N_L_test) = histc(d,bins);
                    
                    
                    d2 = conv(X1_0,ss_n);
                    d2 = d2(T_range_downsampled(1):T_range_downsampled(2));
                    s_all_hh(:,nn,kk,N_L_test)  = histc(d2,bins);

                end
            
                %%
                for N_L_fit = 1:(N_L+1)
                    fprintf('      N_L_fit: %d / %d\n',N_L_fit,N_L+1);
                    %%
                    parfor hh = 1:N_H
                        fprintf('        hh: %d / %d\n',hh,N_H);

                        %%
                        kk_n = stimBasis*glm_k_stim(:,nn,kk,N_L_fit,hh);
                        kk_n = kk_n./sqrt(sum(kk_n.^2));
                        hh_n = spkHistBases_c{hh}*glm_h_spk(1:size(spkHistBases_c{hh},2),nn,kk,N_L_fit,hh);

                        sts_c = sts_glm{N_L_test,N_L_fit,hh}(:);
                        sts_c = sts_c(sts_c >= T_range_downsampled(1) & sts_c <= T_range_downsampled(2));


                        dc = diff(sts_c);
                        if(~isempty(dc))
                            glmStats_isis(:,N_L_test,nn,kk,N_L_fit,hh) =  histc(dc,isi_bins);
                        else
                            glmStats_isis(:,N_L_test,nn,kk,N_L_fit,hh) = 0;
                        end

                        if(~isempty(sts_c))

                            N_spk = length(sts_c);
                            r_bar_glm(nn,kk,N_L_test,N_L_fit,hh) = N_spk./(T_range_downsampled(2)-T_range_downsampled(1) + 1)*1e3;

                            X_spk = X1_0(sts_c + DT_STA + (0:-1:(-T_STA+1)));

                            ss = mean(X_spk,1)';
                            sta_glm(:,nn,kk,N_L_test,N_L_fit,hh) = ss;

                            ss_n = ss./sqrt(sum(ss.^2));





                            sta_n_glm(:,nn,kk,N_L_test,N_L_fit,hh) = ss_n;

                            d = X_spk*ss_n;
                            s_spk_glm(: ,nn,kk,N_L_test,N_L_fit,hh) = histc(d,bins);


                            d2 = conv(X1_0,ss_n);
                            d2 = d2(T_range_downsampled(1):T_range_downsampled(2));
                            s_all_glm(:,nn,kk,N_L_test,N_L_fit,hh)  = histc(d2,bins);
                        end
                    end
                end
            end
        end

    end



    %%
    fits_gain.glm.s_spk = s_spk_glm;
    fits_gain.glm.s_all = s_all_glm;
    fits_gain.glm.sta_n = sta_n_glm;
    fits_gain.glm.sta   = sta_glm;
    fits_gain.glm.r_bar   = r_bar_glm;
    fits_gain.glm.p_spk   = s_spk_glm./sum(s_spk_glm,1);
    fits_gain.hh.s_spk  = s_spk_hh;
    fits_gain.hh.s_all  = s_all_hh;
    fits_gain.hh.sta_n  = sta_n_hh;
    fits_gain.hh.sta   = sta_hh;
    fits_gain.hh.r_bar   = r_bar_hh;
    fits_gain.hh.p_spk   = s_spk_hh./sum(s_spk_hh,1);
    fits_gain.bins = bins;
    
    fits_gain.glm.kh_all = kh_all_glm;
    fits_gain.glm.kh_spk = kh_spk_glm;



    fits_gain.hh.isis = hh_isis;
    fits_gain.glm.isis = glmStats_isis;
    fits_gain.isi_bins = isi_bins;

    fits_gain.alpha = alpha;
    %%

    fprintf('saving...\n');
    save(sprintf('%s/Mease_gainInfo.mat',saveDir),'fits_gain','-v7.3');
    fprintf('done.\n');


end

