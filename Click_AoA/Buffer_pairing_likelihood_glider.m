function L=Buffer_pairing_likelihood_glider(Features,Buffer_Params)

    %% input Features: [ICI, Orientation, Spatial, DTW, Spectral]
    constraint_flag=0;
    % L_prev=0;
    Co=[0.2 0.2 0.6];
    for j=1:3
        x=Features(j);
        if j==2
            objA=Buffer_Params.Orientation.object;
            N_val=Buffer_Params.Orientation.normalize;
        elseif j==1
            objA=Buffer_Params.ICI.object;
            N_val=Buffer_Params.ICI.normalize;
        elseif j==3
           mu = 0;        % mean
           sigma = 5; %20;     % standard deviation
           objA = makedist('Normal', 'mu', mu, 'sigma', sigma);
           N_val=0.0199; 
            constraint_flag=x>60; 
        end
        GMM(j)=Co(j)*pdf(objA,x)/N_val;
    end
    L=sum(GMM);
    if constraint_flag
        L=eps;
    end
end


