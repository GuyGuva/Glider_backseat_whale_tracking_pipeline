function  Whale_output=AoA_Extract(Audio_name,Raw_audio,Click_detections,P,Fs,Plot_flag)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           November 2025
%   DESCRIPTION:

%   This function gets a buffer from a raw audio file with the location 
%   of click detections and attributes them with the corresponding whale source
%   and estimates their angle of arrival with respect to global reference. 

%   INPUT:
%   > Audio_name               - String noting the Audio name
%   > Raw_audio                - Matrix of MX4 comprising 4 channels of M samples of a segment from a the audio file
%   > Click_Detections         - Struct containing the arrival times and amplitudes of detected clicks 
%   > Fs                       - Scalar representing the sampling frequency.
%   > P                        - Struct containing the detector parameters
%   > Plot_flag                - Flag for visualizing results

%   OUTPUT:
%   > Whale_output  - Struct containing the whales and glider parameters at each click arriaval time: 

%   > Whale_output.click_arrival_times_datetime
%   > Whale_output.click_arrival_times_sec
%   > Whale_output.click_amplitudes
%   > Whale_output.glider_latitude
%   > Whale_output.glider_longitude
%   > Whale_output.glider_heading
%   > Whale_output.whale_bearing_from_glider
%   > Whale_output.Glider_depth
%   > Whale_output.Rank (clicks average amplitide X number of clicks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Prepare audio data for analysis
    test_multi=Raw_audio(:,P.Channels_idx);
    test_multi = filtfilt(P.sos,P.g,test_multi);
    t_test=(0:1/Fs: (1/Fs)*(length(test_multi)-1))';
    [start_time,record_time]=Extract_datetime(Audio_name,t_test);
    locs_all=Click_detections.ToAs;

%% run click separation algorithm 
    t_index=0;
    NOB=round(t_test(end)/P.Buffer_length); % Number of buffers
    for Buf_ind=1:NOB
         disp(['Separating clicks: ' num2str(100*Buf_ind/NOB) '%'])
        t1=(Buf_ind-1)*P.Buffer_length;
        if Buf_ind<NOB
            Buffer_lims=(Buf_ind-1)*P.Buffer_length*Fs+1:Buf_ind*P.Buffer_length*Fs;
        else
            Buffer_lims=(Buf_ind-1)*P.Buffer_length*Fs+1:t_test(end)*Fs;
        end
        locs=locs_all(locs_all>t1 & locs_all<t1+P.Buffer_length);
        whale_clicks=Extract_DoAs_per_buffer(locs,start_time,t_test,P,t1,test_multi,Buffer_lims,Fs,Plot_flag);
        Whales_buffer(Buf_ind).whale_clicks=whale_clicks;
        t_index=t_index+1;
        for i=1:length(Whales_buffer(Buf_ind).whale_clicks)
            Detections(t_index).ToAs(i)={Whales_buffer(Buf_ind).whale_clicks(i).ToAs'};
            Detections(t_index).Pkk(i)={Whales_buffer(Buf_ind).whale_clicks(i).Pks};
            Az=Whales_buffer(Buf_ind).whale_clicks(i).Azimuth_to_whale;
            El=Whales_buffer(Buf_ind).whale_clicks(i).Elevation_to_whale;
            if isempty(Az) | isempty(El)
                Detections(t_index).DoA(i)={1e3*rand(1,2)};
                Detections(t_index).DoA_all_clicks(i)={1e3*rand(1,2)};
            else
                Detections(t_index).DoA(i)={[mean(Az) mean(El)]};
                Detections(t_index).DoA_all_clicks(i)={[Az' El']};
            end
            Detections(t_index).ICI(i)={diff(sort(Whales_buffer(Buf_ind).whale_clicks(i).ToAs'))};
        end
    
    end
    t_index=t_index+1;
    Detections(t_index).ToAs={};
    Detections(t_index).Pkk={};
    Detections(t_index).DoA={};
    Detections(t_index).DoA_all_clicks={};
    Detections(t_index).ICI={};
  
%%  Click train formation (sequences' association between buffers)
    [trajectories,id_j_ToAs,DoA_trajectories]=run_train_DoA_version(Detections,test_multi(:,1),Fs,P,Plot_flag);

%% Association of click trains
    Whale_output=Associate_trains_glider_V2(P,record_time,DoA_trajectories,trajectories,id_j_ToAs,test_multi(:,1),Fs,Plot_flag);

end