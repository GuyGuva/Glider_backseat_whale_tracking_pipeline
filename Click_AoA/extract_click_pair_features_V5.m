function [alpha_AoA,alpha_Pkk,DTW,Membership]=extract_click_pair_features_V5(F_weights,ICI,c1,c2,AoA_vec1,AoA_vec2,F_ds,All_objs)
  
    % c1=Click_i; c2=Click_j; AoA_vec1=AoA_all(i,:); AoA_vec2=AoA_all(j,:);


        Pkk_left=max(c1)-min(c1);
        Pkk_right=max(c1)-min(c2);
       

        DTW=dtw(c1/max(c1),c2/max(c2));
        alpha_Pkk=log(Pkk_right/Pkk_left);
        Spectral=spectral_similarity(c1,c2,F_ds,0);
        alpha_AoA=norm(AoA_vec2-AoA_vec1);
        Features=[alpha_AoA alpha_Pkk DTW Spectral];
        sigma = 4;
        mu = 0;              

        for q=1:4
            if q==1
                Mem_val(q)=(1/(0.08*sigma*sqrt(2*pi))) * exp(-0.5 * ((alpha_AoA - mu)/sigma).^2);
            else
                objA=All_objs(q).objs;
                N_val=All_objs(q).N_val;
                Mem_val(q)=pdf(objA,Features(q))/N_val;
            end
        end
       % Membership=sum(Mem_val);
       % F_weights=[1 0 0 0];
       % F_weights=[0.7 0.1 0.1 0.1];
       Membership=sum(F_weights.*Mem_val);
       
end

