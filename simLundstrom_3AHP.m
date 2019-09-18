function [spikeTimes,V] = simLundstrom_3AHP(i_stim,h,i_mult,i_dc,stim_gain)

    if(nargin < 5)
        stim_gain = @(t) ones(size(t));
    end
    
    T_int = length(i_stim);

    if(length(i_mult) ~= length(i_dc))
        error("i_mult and i_dc must be the same length.");
    end
    
    N = length(i_mult);
    if(N==1)
        N = length(i_dc);
    end

    if(N==1)
        N = length(stim_gain(1));
    end

    minSpikeInterval = 2.0;
    spikeThreshold   = -10.0;

    lastSpike =  zeros(N,1);


    nParams = 7;
    returnV = nargout > 1;
    V_c = zeros(nParams,N);
    V_c(1,:) = -70.0;
    V_c(2,:) = m_inf(V_c(1,:));
    V_c(3,:) = h_inf(V_c(1,:));
    V_c(4,:) = n_inf(V_c(1,:));
    if(returnV)
        V = zeros(nParams,N,T_int);
        V(:,:,1) = V_c;
    end



    spikeTimes = cell(N,1);
    for jj = 1:N
        spikeTimes{jj} = [];
    end



    dispInterval = max(1,floor(T_int*0.01));
    fprintf("Simulation beginning...\n");
    
%     k1 = zeros(nParams,N);
%     k2 = zeros(nParams,N);
%     k3 = zeros(nParams,N);
%     k4 = zeros(nParams,N);
    
    stim = i_stim(1)*stim_gain(1).*i_mult+i_dc;
    previouslyCrossed = false(1,N);
    for t_idx = 1:(T_int-1)

        if(mod(t_idx , dispInterval)==0)
            fprintf("Simulation %d%% complete.\n",round(((t_idx*1.0)/T_int)*100.0));
        end
        
        
        k1 = computeHHvals(V_c,stim);

        k2 = computeHHvals(V_c + h/2.0*k1,stim);

        k3 = computeHHvals(V_c + h/2.0*k2,stim);

        stim = i_stim(t_idx+1)*stim_gain(t_idx+1).*i_mult+i_dc;
        k4 = computeHHvals(V_c + h*k3,stim);

        V_c = V_c + (k1 + 2.0.*k2 + 2.0.*k3 + k4).*(h./6.0);

        V_c(1,:) = max(min(V_c(1,:),150.0),-250.0);
        V_c(2,:)  = max(min(V_c(2,:),1.0),0.0);
        V_c(3,:)  = max(min(V_c(3,:),1.0),0.0);
        V_c(4,:)  = max(min(V_c(4,:),1.0),0.0);
        V_c(5,:)  = max(V_c(5,:),0.0);
        V_c(6,:)  = max(V_c(6,:),0.0);
        V_c(7,:)  = max(V_c(7,:),0.0);
        
        currentlyCrossed = (V_c(1,:) >= spikeThreshold);

        for jj = 1:N
            if(currentlyCrossed(jj) && ~previouslyCrossed(jj))
                if(((t_idx-lastSpike(jj))*h >= minSpikeInterval) || lastSpike(jj) == 0)
                    lastSpike(jj)   = t_idx;
                    spikeTimes{jj} = [spikeTimes{jj};t_idx];

                    V_c(5:7,jj) = V_c(5:7,jj) + 1;
                end
            end
        end
        if(returnV)
            V(:,:,t_idx+1) = V_c;
        end
        previouslyCrossed = currentlyCrossed;
    end
    fprintf("Simulation complete.\n");

end


function [m] = m_inf(v)
    m = 0.1*(v+40) ./( 0.1*(v+40)  + 4*exp(-(v+65)./18.0).*(1-exp(-0.1*(v+40))));
end

function [n] = n_inf( v)
    n =  0.01*(v+55)./(0.01*(v+55) + 0.125*exp(-(v+65)./80).*(1-exp(-0.1*(v+55))));
end

function [h] = h_inf(v)
    h = 0.07./(0.07 +exp((v+65)./20)./(1+exp(-0.1*(v+35))));
end

function [x_c] = computeHHvals(x_c,i_stim)

    E_l = -54.4;
    E_Na = 50.0;
    E_K  = -77.0;
    E_r  = -100;

    V = x_c(1,:);
    m = x_c(2,:);
    h = x_c(3,:);
    n = x_c(4,:);
    a1 = x_c(5,:);
    a2 = x_c(6,:);
    a3 = x_c(7,:);


    g_bar_na = 120;
    g_bar_k = 36;

    %CAP = 1.0;

    G_Na = g_bar_na*max(0,min(1,m.^3 .* h));
    G_K = g_bar_k*max(0,min(1,n.^4));
    G_l = 0.3;
    G_a = G_l*max(0,a1*0.05 + a2*0.006 + a3*0.004);

    x_c(1,:)=(G_l*E_l+i_stim(:)' + G_Na*E_Na + G_K*E_K + G_a*E_r  -V.*(G_l + G_Na + G_K+G_a));
    
    x_c(2,:)=((0.1*(V+40) ./( 0.1*(V+40)  + 4*exp(-(V+65)./18.0).*(1-exp(-0.1*(V+40))))) - m)./( 1.0./( 0.1*(V+40)./(1-exp(-0.1*(V+40)))  + 4*exp(-(V+65)./18.0)));
    x_c(3,:)=((0.07./(0.07 +exp((V+65)./20)./(1+exp(-0.1*(V+35))))) - h)./(1.0./(0.07*exp(-(V+65)./20) +1./(1+exp(-0.1*(V+35)))));
    x_c(4,:)=((0.01*(V+55)./(0.01*(V+55) + 0.125*exp(-(V+65)./80).*(1-exp(-0.1*(V+55))))) - n)./(1.0./(0.01*(V+55)./(1-exp(-0.1*(V+55))) + 0.125*exp(-(V+65)./80)));
    

    
    x_c(5,:)=-a1./300.0;
    x_c(6,:)=-a2./1000.0;
    x_c(7,:)=-a3./6000.0;
end





