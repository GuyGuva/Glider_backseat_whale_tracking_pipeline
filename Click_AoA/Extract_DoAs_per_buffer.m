function  whale_clicks=Extract_DoAs_per_buffer(locs,start_time,t_test,P,t1,test_multi,Buffer_lims,Fs,Plot_flag)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           November 2025
%   DESCRIPTION:

%   This function gets a short buffer (3 sec) from a filtered audio file with the location 
%   of click detections and attributes them with the corresponding whale source
%   and estimates their angle of arrival with respect to global reference. 

%   INPUT:
%   > test_multi               - Matrix of MX4 comprising 4 channels of M samples of a segment from a the audio file
%   > start_time               - Vector of MX1 representing of M datetime samples of a segment from a the audio file
%   > t_test                   - Vector of MX1 representing of M time samples [0,M/Fs] of a segment from a the audio file
%   > t1                       - Scalar indicating the start time of the buffer with respect to the entire recording
%   > Buffer_lims              - Vector of KX1 (k=P.Buffer_length*Fs) representing k time samples of a segment from a the audio file
%   > locs                     - Vector of 1XN with time of arrival [in seconds] for N most dominant transients in the buffer
%   > Fs                       - Scalar representing the sampling frequency.
%   > P                        - Struct containing the detector parameters
%   > Plot_flag                - Flag for visualizing results

%   OUTPUT:
%   > Whale_output  - Struct containing the whales and glider parameters at each click arriaval time: 

%   > Whale_output.ToAs
%   > Whale_output.Pks
%   > Whale_output.Azimuth_to_whale
%   > Whale_output.Elevation_to_whale
%   > Whale_output.Glider_heading
%   > Whale_output.Glider_Lat
%   > Whale_output.Glider_Lon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
    whale_clicks.ToAs=[];
    whale_clicks.Pks=[];
    whale_clicks.Azimuth_to_whale=[];
    whale_clicks.Elevation_to_whale=[];
    whale_clicks.Glider_heading=[];
    whale_clicks.Glider_Lat=[];
    whale_clicks.Glider_Lon=[];
    if length(locs)>P.NOT
    %% Estimate the AoA for each detected click:
        Locs=locs-t1;
        Locs(Locs<P.W)=[];
        Locs(Locs>P.Buffer_length-P.W)=[];
    
        buffer_all_channels=test_multi(Buffer_lims,:);
        Locs_datetime=start_time(int32((Locs+t1)'*Fs));
        Orientation=Extract_glider_orientation(P.AllData,Locs_datetime);
        locs=Locs;
    
        bearing_correct=zeros(1,length(locs));
        elevation_correct=zeros(1,length(locs));
        theta_body_est=zeros(1,length(locs));
        phi_body_est=zeros(1,length(locs));
        for click_idx=1:length(locs)
           Glider_sensor_data.heading=Orientation.heading_at_clicks(click_idx);
           Glider_sensor_data.pitch=Orientation.pitch_at_clicks(click_idx);
           Glider_sensor_data.roll=Orientation.roll_at_clicks(click_idx);
           Glider_sensor_data.declination=Orientation.declination_at_clicks(click_idx);
        
           Click_ToA=Locs_datetime(click_idx);
           Click_ToA.TimeZone = 'UTC';
    
           ROI= buffer_all_channels(Fs*(locs(click_idx)-P.roi):Fs*(locs(click_idx)+P.roi),:)';
    
           [bearing_correct(click_idx), elevation_correct(click_idx),theta_body_est(click_idx), phi_body_est(click_idx)]=Estimate_angles(P.sensor_pos,Glider_sensor_data,ROI,Fs);
    
        end
        AoA_all=[bearing_correct' elevation_correct'];
            
    %% Cluster clicks
        Nd=length(locs);
        L_tot = 1e9*ones(Nd,Nd);
        ICI=zeros(Nd,Nd);
        for i=1:Nd-1
               Click_i= buffer_all_channels(Fs*(locs(i)-P.W):Fs*(locs(i)+P.W),P.channel);
           for j=i+1:Nd
               Click_j= buffer_all_channels(Fs*(locs(j)-P.W):Fs*(locs(j)+P.W),P.channel);
                if locs(j)-locs(i)>P.ICI_min && locs(j)-locs(i)<P.ICI_max 
                   ICI(i,j)= locs(j)-locs(i);
                    [~,~,~,Membership]=extract_click_pair_features_V5(P.F_weights,ICI(i,j),Click_i,Click_j,AoA_all(i,:),AoA_all(j,:),Fs,P.All_objs);                   
                    L_tot(i,j)=-log(Membership);                     
                end          
           end
        end
        [C,~] = do_assignment (L_tot,-log(-P.lone_p_click_separaton), 1e9);
        nchain = 0;
        chain_nr = zeros(Nd,1);
        for i = 1:Nd
            j = find(C(i,1:Nd));
            if (~isempty(j))
                if (chain_nr(i) && ~any(chain_nr(j)))
                    chain_nr(j) = chain_nr(i);
                else
                    if(chain_nr(i)==0 && ~any(chain_nr(j)))
                        nchain = nchain + 1;
                        chain_nr(j) = nchain;
                        chain_nr(i) = nchain;
                    end
                end
            end
        end
        All_traces=chain_nr;   
        All_traces(All_traces==0)=[];
        Detected_traces=unique(All_traces)';
         Gr={};
         if ~isempty(Detected_traces)
            for i=Detected_traces
                    Gr(i)={find(chain_nr==i)};
            end            
            Detected_subtrains= Gr(~cellfun(@isempty, Gr));
         else
             Detected_subtrains={};
         end
         Detected_subtrains = Detected_subtrains(cellfun(@(x) numel(x) >= 3, Detected_subtrains));
    
     %% Remove click groups with unstable AoA
        Azi=bearing_correct;
        Ele=elevation_correct;
        Var=zeros(1,length(Detected_subtrains));
        for i=1:length(Detected_subtrains)
             subtrain_idx=Detected_subtrains{i};
             [~,TFrm] = rmoutliers(Azi(subtrain_idx)); Azi(subtrain_idx(TFrm))=median(Azi(subtrain_idx));
             [~,TFrm] = rmoutliers(Ele(subtrain_idx)); Ele(subtrain_idx(TFrm))=median(Ele(subtrain_idx));       
             Var(i)=mean(std([Azi(subtrain_idx)' Ele(subtrain_idx)']));
        end
        Include=Var<P.angle_variance;
        final_subtrains=Detected_subtrains(Include);
    
     %% Visualize results 
        co={'r*','g*','c*','k*','y*','m*','b*','rx','bx','cx','kx','mx','gx','bx','ro','go','co','ko','mo','yo','bo'};
        if Plot_flag.Click_AoA_per_buffer
            figure;
            subplot(1,3,1);
            plot(t_test,test_multi(:,1)); hold on;
            for i=1:length(final_subtrains)
                ind=final_subtrains{i};     
                Pk=Peaks_extract(test_multi(:,1),t1+locs,Fs);
                dotH=plot(t1+locs(ind),Pk(ind),[co{i}],'LineWidth',2); hold on; 
                   pause(0.05);  % calls DRAWNOW implicitly
                   set(dotH, 'XData', t1+locs(ind), 'YData', Pk(ind)); hold on;
            end
            % xlim([t1 t1+3]);
        end   
    
        if Plot_flag.Click_AoA_per_buffer
           subplot(1,3,2);
        end
        Bearings_to_whale=[];
        for i=1:length(final_subtrains)
             subtrain_idx=final_subtrains{i};
             Bearings_to_whale=[Bearings_to_whale mean(Azi(final_subtrains{i}))];
             if Plot_flag.Click_AoA_per_buffer
                 plot(Azi(subtrain_idx),Ele(subtrain_idx),[co{i}],'Linewidth',2); grid on; hold on;
             end
        end
        if Plot_flag.Click_AoA_per_buffer
            xlabel('Azimuth [\circ]'); ylabel('Elevation [\circ]');
            axis([-180 180 -180 180]);
        end  
        if Plot_flag.Click_AoA_per_buffer
            lat = Orientation.Lat_at_clicks(1);
            lon = Orientation.Lon_at_clicks(1); 
            visualize_on_map(1e3,Bearings_to_whale,co,lat,lon,mean(Orientation.heading_at_clicks));
        end
    
       %% Store results
        for i=1:length(final_subtrains)
            ind=final_subtrains{i}; 
            Pk=Peaks_extract(test_multi(:,1),t1+locs,Fs);
            whale_clicks(i).ToAs=t1+locs(ind);
            whale_clicks(i).Pks=Pk(ind);
            whale_clicks(i).Azimuth_to_whale=Azi(ind);
            whale_clicks(i).Elevation_to_whale=Ele(ind);
            whale_clicks(i).Glider_heading=Orientation.heading_at_clicks(ind);
            whale_clicks(i).Glider_Lat=Orientation.Lat_at_clicks(ind);
            whale_clicks(i).Glider_Lon=Orientation.Lon_at_clicks(ind);
        end
    end
           
end

