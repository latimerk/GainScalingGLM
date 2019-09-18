function [cycle] = cycleAverageSpikes(NT_period,sts,T_range,StimPeriods)

    %NT_period = 40;

    NP = size(sts,1);
    NL = size(sts,2);
    cycle = zeros(NT_period,NP,NL);

    for ii = 1:NP
        TL = StimPeriods(ii)*1000;

        ps = T_range(1):TL:T_range(1);

        NC = length(ps)-1;

        bins = linspace(0,TL,NT_period+1);

        for jj = 1:NL
            for kk = 1:NC
                for ll = 1:NT_period
                    T_start = ps(kk) + bins(ll);
                    T_end   = ps(kk) + bins(ll+1);
                    cycle(ll,ii,jj) = cycle(ll,ii,jj) + sum((sts{ii,jj} >= T_start) & (sts{ii,jj} < T_end))./(1e-3*(T_end-T_start));
                end
            end
        end
        cycle(:,ii,:) = cycle(:,ii,:)./NC;

    end
end