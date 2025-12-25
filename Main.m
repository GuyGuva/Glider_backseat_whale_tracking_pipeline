function Whale_output=Main(File_directory,File_name)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Define directories path
     Main_folder=pwd; 
     Audio_folder=File_directory; 
     Program_folder.Click_Detector=[Main_folder '/Click_Detector'];
     Program_folder.Click_AoA=[Main_folder '/Click_AoA'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Parameter settings
    
     % User defined parameters:
     Buffer_size=30;             % Set the length of the analyzed window (in [sec])
     channel_select=1;           % Select channel for detection analysis
     Buffer_index=0;             % Defines the buffer number within the audio
     D_Threshold=1.42;           % set the click presence detection threshold
     SNR_thresh=46;              % Min SNR threshold [in dB] for remote recordings
     Whale_output=[];

     % Default parameters:
     cd(Program_folder.Click_Detector);
     Presence_detector_settings=Default_parameters_presence_detector(SNR_thresh,D_Threshold);
     cd(Program_folder.Click_AoA);
     Click_Separator_settings=Default_parameters_click_separator;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Visualization
      % Plot flags for testing functions
         Plot_flag.Click_presence_detection=0;
         Plot_flag.Click_detection=0;
         Plot_flag.Click_AoA_per_buffer=0;
         Plot_flag.Click_AoA_trains_formation=0;
         Plot_flag.Click_AoA_final=1;
         Plot_flag.Click_AoA_final_time_series=0;
         Plot_flag.Click_AoA_final_time_series_smoothed=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %% Select audio file
    Audio_name=File_name;    
    while(1) 
        Buffer_index=Buffer_index+1;
        cd(Program_folder.Click_Detector);
    
        %% Block 1: obtain buffer
        [Raw_audio,Fs,Edge_flag]=read_audio(Audio_folder,Program_folder.Click_Detector,Audio_name,Buffer_size,Buffer_index);  
        if Edge_flag | isempty(Raw_audio)
            break;
        end

        %% Block 2: Click detection    
         [Click_detections,Echolocation_clicks_Presence_flag]=Click_Detector(Raw_audio(:,channel_select)',Presence_detector_settings,Fs,Buffer_size,Plot_flag);
    
        %% Block 3: Click AoA estimation
        if Echolocation_clicks_Presence_flag
            disp('Clicks are found');
            cd(Program_folder.Click_AoA);
            Whale_output=AoA_Extract(Audio_name,Raw_audio,Click_detections,Click_Separator_settings,Fs,Plot_flag);
        else
            disp('No clicks are found');
        end
    end

    cd(Main_folder)

end