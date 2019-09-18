function [bestStim,rs,bestStim_interp] = getStimulusRate(searchType,searchRange,targetRate,sigMuRatio,T_range,Stim,h,G_Na,G_K)
    rs = nan(length(searchRange),1);
    fprintf("Searching for target rate.\n");

    if(strcmpi(searchType,"Lundstrom3AHP"))
        S = simLundstrom_3AHP(Stim,h,sigMuRatio*searchRange,searchRange);
    elseif(strcmpi(searchType,"Lundstrom2Na"))
        S = simLundstrom_2Na(Stim,h,sigMuRatio*searchRange,searchRange);
    elseif(strcmpi(searchType,"LundstromDefault"))
        S = simLundstrom_Default(Stim,h,sigMuRatio*searchRange,searchRange);
    elseif(strcmpi(searchType,"Mease"))
        S = simMease(Stim,h,sigMuRatio*searchRange,searchRange,G_Na,G_K);
    end
    for ii = 1:length(searchRange)
        if(isempty(S{ii}))
            rs(ii) = 0;
        else
            nspks = sum((S{ii}>=T_range(1)) & (S{ii} <= T_range(end)));
            rs(ii) = nspks/((T_range(end) - T_range(1) + 1.0)*(h/1000.0));
        end

    end
    [~,idx] = min(abs(rs - targetRate));
    bestStim = searchRange(idx);
    

    r1 = find(rs <= targetRate,1,'last');
    r2 = find(rs >= targetRate,1,'first');
    if(~isempty(r1) && ~isempty(r2))
        if(r2 == r1+1)
            dy = rs(r2) - rs(r1);
            fr = (targetRate-rs(r1))/dy;
            
            dx = searchRange(r2) - searchRange(r1);
            
            bestStim_interp = dx*fr + searchRange(r1);
            
        else
            bestStim_interp = bestStim;
        end
    else
        bestStim_interp = bestStim;
    end
    
    fprintf(" found stim level %.2f giving rate %.2f (interp rate = %.3f)\n", bestStim, rs(idx), bestStim_interp);

end