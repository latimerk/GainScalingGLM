function [spikeTimes,V,S] = simMease(i_stim,h,i_mult,i_dc,G_Na,G_K,stim_gain)

    if(nargin < 7)
        stim_gain = @(t) ones(size(t));
    end
    
    T_int = length(i_stim);

    
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


    nParams = 4;
    returnV = nargout > 1;
    returnS = nargout > 2;
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
    
    
    t_idx = 1;
    %i_stim = i_stim.*stim_gain(1:T_int)';
    stim = i_stim(t_idx,:)*i_mult+i_dc;
    previouslyCrossed = false(1,N);
    
    
    if(returnS)
        S = zeros(N,T_int);
        S(:,1) = stim;
    end
    
    Nspks = zeros(N,1);

    for t_idx = 1:(T_int-1)

        if(mod(t_idx , dispInterval)==0)
            idx = 1;%find(i_dc == 0.35,1);
            fprintf("Simulation %d%% complete. (rate %d = %.2f)\n",round(((t_idx*1.0)/T_int)*100.0), idx,Nspks(idx)/(t_idx*h/1e3));
        end
        
        
        k1 = computeHHvals(V_c,stim,G_Na,G_K);

        k2 = computeHHvals(V_c + h/2.0*k1,stim,G_Na,G_K);

        k3 = computeHHvals(V_c + h/2.0*k2,stim,G_Na,G_K);

        stim = i_stim(t_idx+1,:)*i_mult+i_dc;
        k4 = computeHHvals(V_c + h*k3,stim,G_Na,G_K);

        V_c = V_c + (k1 + 2.0.*k2 + 2.0.*k3 + k4).*(h./6.0);

        V_c(1,:) = max(min(V_c(1,:),150.0),-250.0);
        V_c(2,:)  = max(min(V_c(2,:),1.0),0.0);
        V_c(3,:)  = max(min(V_c(3,:),1.0),0.0);
        V_c(4,:)  = max(min(V_c(4,:),1.0),0.0);
        

        currentlyCrossed = (V_c(1,:) >= spikeThreshold);
        for jj = 1:N
            if(currentlyCrossed(jj) && ~previouslyCrossed(jj))
                if(((t_idx-lastSpike(jj))*h >= minSpikeInterval) || lastSpike(jj) == 0)
                    lastSpike(jj)   = t_idx;
                    spikeTimes{jj} = [spikeTimes{jj};t_idx];
                    Nspks(jj) = Nspks(jj) + 1;
                end
            end
        end
        if(returnV)
            V(:,:,t_idx+1) = V_c;
        end
        if(returnS)
            S(:,t_idx+1) = stim;
        end
        previouslyCrossed = currentlyCrossed;
    end
    fprintf("Simulation complete.\n");

end

function aa = A_func(A, K, v, th)
   aa = A.*K*ones(size(v));
   for ii = 1:length(v)
       if(abs(v(ii)-th) > 1e-8)
           aa(ii) = (A.*(v(ii)-th)./(1.0-exp(-(v(ii)-th)./(K))));
       end
   end
end

function B = B_func( A, K, v,  th)
    B =  A_func(A,K,-v,-th);
end

function C = C_func(V, A1, A2, K, th)
    C = 1.0./(1.0+exp(-(V-th-(K.*log(A2./A1)))./K));
end

function m = m_inf(V)
    m = C_func(V,0.182,0.124,9.0,-35.0);
end
function h = h_inf(V)
    h = 1.0./(1.0+exp((V+65.0)./6.2));
end
function n = n_inf(V)
    n = C_func(V,0.020,0.002,9.0,20.0);
end

function t = taui_m(V,tadj)
    t = ((A_func(0.182,9.0,V,-35.0)+B_func(0.124,9.0,V,-35.0))*tadj);
end
function t = taui_n( V,tadj)
    t = ((A_func(0.020,9.0,V,20.0)+B_func(0.002,9.0,V,20.0))*tadj);
end
function t = taui_h(V,tadj)
    t = ((A_func(0.024 ,5.0,V,-50.0)+B_func(0.0091,5.0,V,-75.0))*tadj);
end


function [x_c] = computeHHvals(x_c,i_stim_c,g_bar_na,g_bar_k)

    E_l = -70.0;
    E_Na = 50.0;
    E_K  = -77.0;

    V = x_c(1,:);%min(max(,-250.0),150.0)
    m = x_c(2,:);
    h = x_c(3,:);
    n = x_c(4,:);

    %celsius = 25;
    
    t_adj_n  = 1;%2.3^((celsius - 16)/10);%3.2094;
    t_adj_m  = 1;%2.3^((celsius - 23)/10);
    g_bar_na = 0.1*g_bar_na;%t_adj*0.1*g_bar_na;
    g_bar_k  = 0.1*g_bar_k;%t_adj*0.1*g_bar_k;

    %CAP = 1.0;

    G_Na = g_bar_na*max(0,min(1,m.^3 .* h));
    G_K = g_bar_k*max(0,min(1,n));
    G_l = 0.04; %25 ms leak! Paper said 25ms but has 0.025 conductance

%     x_c(1,:)=((G_l*(E_l -V) + G_Na.*(E_Na - V) + G_K .*(E_K -V) + G_a.*(E_r-V) + i_stim(:)'));%./CAP
    x_c(1,:)=(G_l*E_l+i_stim_c(:)' + G_Na*E_Na + G_K*E_K   -V.*(G_l + G_Na + G_K));%./CAP
    
    x_c(2,:)= (m_inf(V) - m).*taui_m(V,t_adj_m);
    x_c(3,:)= (h_inf(V) - h).*taui_h(V,t_adj_m);
    x_c(4,:)= (n_inf(V) - n).*taui_n(V,t_adj_n);
    

%     x_c(2,:)=(1.0./(1.0+exp(-(V+35.0-(-3.453526))./9.0))  - m).*((A_func(0.182,9.0, V,-35.0)+A_func(0.124,9.0, -V, 35.0))*3.2094);
%     x_c(3,:)=(1.0./(1.0+exp((V+65.0)./6.2))               - h).*((A_func(0.024 ,5.0,V,-50.0)+A_func(0.0091,5.0,-V, 75.0))*3.2094);
%     x_c(4,:)=(1.0./(1.0+exp(-(V-20.0-(20.723266))./9.0))  - n).*((A_func(0.020,9.0, V, 20.0)+A_func(0.002,9.0, -V,-20.0))*3.2094);


end





