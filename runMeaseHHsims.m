%Simulate HH neurons as in
% Mease, R. A., Famulare, M., Gjorgjieva, J., Moody, W. J., and Fairhall, A. L. (2013). Emergence of adaptive computation by single neurons in the developing cortex.
%   Journal of Neuroscience
%

rng(08012018);

debug = false;






if(debug)
    BUFFER = 1;
    L = 1000*100-2;
    L_init = L;
else
    BUFFER = 60*1000*100;
    L = 2000*1000*100;
    L_init = 100*1000*100;
end
downsampleRate = 100;
h = 1.0/(downsampleRate);
T_range = BUFFER+[1 L];
T_range_init = BUFFER+[1 L_init];
targetRate = 10.0;
maxZeroRate= 5;

xs = -0.1:0.01:3;
sigMuRatio = 4.0;
StimLevels = [1.0 1.3 1.6 2.0];

frFunctions = nan(length(xs),15,15);


T_range_downsampled = ceil(T_range./downsampleRate);



X1 = randn(ceil((L+2*BUFFER)/downsampleRate),1);

y = zeros(length(X1)*downsampleRate,1);
for ii = 1:length(X1)
    y((1:downsampleRate) + (ii-1)*downsampleRate) = X1(ii);
end

y_init = y(1:(T_range_init(end) + 1));

G_mat_info_paper = [1 5;
              1 7;
              1 8;
              1 10;
              1 11;
              2 12;
              3 14;
              3 15;
              4 15;
              5 15;
              5 15;
              6 15;
              7 15;
              8 15;
              9 15];

G_mat_info = G_mat_info_paper;
G_mat_info(:,1) = 1;
G_mat_info(:,end) = 15;

NC = size(G_mat_info,1);
Gs = 500.0+(1:NC).*100.0;

stim_mu = nan(NC,NC);

sts = cell(NC,NC,length(StimLevels));
sts_0 = cell(NC,NC,length(StimLevels));

for ii = 1:size(sts_0,1)
    for jj = 1:size(sts_0,2)
        for kk = 1:size(sts_0,3)
            sts_0{ii,jj,kk} = -1;
            sts{ii,jj,kk} = -1;
        end
    end
end

G_Na = Gs(1);

%%

fprintf("Starting HH sims.\n")
if(exist('Results/Mease/MeaseSims_part.mat','file'))
    load('Results/Mease/MeaseSims_part.mat');
    aa_start = aa+1;
else
    aa_start = 1;
end
for aa = aa_start:15
    G_Na = Gs(aa);
    fprintf("G_Na level %d (%d).\n",aa,G_Na);
    bb = G_mat_info(aa,1):G_mat_info(aa,2); 


    G_Ks  = Gs(bb);


    stim_mu_c = zeros(length(G_Ks),1);
    fr =  nan(length(xs),length(G_Ks));
    parfor bb_c = 1:length(G_Ks)
        [~,fr(:,bb_c),stim_mu_c(bb_c)] = getStimulusRate("Mease",xs,targetRate,sigMuRatio,T_range_init,y_init,h,G_Na,G_Ks(bb_c));
    end
    
    stim_mu_c2 = nan(1,15);
    stim_mu_c2(bb) = stim_mu_c;
    
    fr2 = nan(length(xs),15);
    fr2(:,bb) = fr;
    
    frFunctions(:,aa,:) = fr2;
    stim_mu(aa,:) = stim_mu_c2;
    
%     frFunctions(:,aa,bb) = fr;
%     stim_mu(aa,bb) = stim_mu_c;

    %StimStds = stim_mu_c*StimLevels*sigMuRatio;


    fprintf("Running simulation battery...\n")
    sts_c = cell([1,length(G_Ks),length(StimLevels)]);
    for bb_c = 1:length(G_Ks)
        for ii = 1:length(StimLevels)
            sts_c{1,bb_c,ii} = -1;
        end
    end
    parfor bb_c = 1:length(G_Ks)
        if(fr(1,bb_c) < maxZeroRate)
        
            sigs = stim_mu_c(bb_c)*StimLevels*sigMuRatio; 
            sts_c(1,bb_c,:) = simMease(y,h,sigs(:),stim_mu_c(bb_c),G_Na,G_Ks(bb_c)); 
        end
    end

    sts_0(aa,bb,:) = sts_c;
    
    for ii = G_mat_info(aa,1):G_mat_info(aa,2)
        for jj = 1:length(StimLevels)
            
            if(~isempty(sts_0{aa,ii,jj}))
                sts{aa,ii,jj} = ceil(sts_0{aa,ii,jj}./downsampleRate);
            else
                sts{aa,ii,jj} = [];
            end
        end
     end
   


    save("Results/Mease/MeaseSims_part.mat","-v7.3","aa","frFunctions","xs","sts","sts_0","stim_mu","sigMuRatio","targetRate","maxZeroRate","downsampleRate","X1","T_range_downsampled","debug","StimLevels");
end

% for aa = 1:15
%      for ii = G_mat_info(aa,1):G_mat_info(aa,2)
%         for jj = 1:length(StimLevels)
%             
%             if(~isempty(sts_0{aa,ii,jj}))
%                 sts{aa,ii,jj} = ceil(sts_0{aa,ii,jj}./downsampleRate);
%             else
%                 sts{aa,ii,jj} = [];
%             end
%         end
%      end
% end


fprintf("Saving...\n")

save("Results/Mease/MeaseSims.mat","-v7.3","aa","frFunctions","xs","sts","sts_0","stim_mu","G_mat_info_paper","sigMuRatio","targetRate","maxZeroRate","downsampleRate","X1","T_range_downsampled","debug","StimLevels");
fprintf("done.\n")
