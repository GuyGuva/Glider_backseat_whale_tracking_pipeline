
function P=Default_parameters_click_separator

%% Parameters settings:
     P.Channels_idx=[2 4 1 3];        % Order of channels with respect of the predefind glider frame
     P.roi=1.5e-3;                    % region of interest (roi)- defines the time window [in sec] around clicks for AoA analysis
     P.Buffer_length=3;               % Analysis buffer length [sec] 
     P.lone_p=5;                      % Lone penalty (used for click train formation)      
     P.ITI_min=3;                     % minimum allowed time gap between trains(Inter-train interval (ITI))  [sec]
     P.ITI_max=30;                    % maximum allowed time gap between trains (Inter-train interval (ITI))  [sec]
     P.ICI_min=0.4;                   % minimum allowed time gap between clicks (Inter-click interval (ICI))  [sec]
     P.ICI_max=1.5;                   % maximum allowed time gap between clicks (Inter-click interval (ICI))  [sec]
     P.angle_variance=5;              % maximum allowed angle variance [in degrees] within a click train (within a buffer) 
     P.max_AoA_change=20;             % maximum allowed angle variance [in degrees] between click traind (between buffers) 
     P.transients_threshold=2;        % minimum number of transient detections for considering further processing
     P.lone_p_click_separaton=-0.1;   % Lone penalty (used for click separation within buffers)
     P.channel=1;                     % Set channel for click separation analysis
     P.min_numer_of_clicks_per_whale=20; % Set minimum accepted numer of clicks per whale
     P.W=3e-3;                        % Set pulse width for attributes extraction
     P.NOT=2;                   % Set maximum number of clicks analyzed for the trains' association step
     P.n_clicks=16;                   % Set minimum number of transient for clustering analysis
     d = 0.15;                        % Set the distance [in m] between sensors in the array
     P.sensor_pos = [
            0,     0,              0;
            d,     0,              0;
            d/2,   d*sqrt(3)/2,    0;
            d/2,   d*sqrt(3)/6,    d*sqrt(2/3)   % Define sensors postion in the array frame
        ];
  %% import trained models
     load All_objs.mat;       % load the GMM parameters of the click sequence similarity attributes (buffer level)
     load Buffer_Params;      % load the GMM parameters of the clicks' similarity attributes (click level)
     load F_weights;          % load the GMM parameters of the clicks' similarity attributes (buffer level)
     load hp_filter_sos.mat;  % load BPF (2-24khz) parameters (sos and g)
     P.sos=sos;
     P.g=g;
     P.Buffer_Params=Buffer_Params;    % - struct containing the GMM parameters of the clicks' similarity attributes in the buffer level
     P.F_weights=F_weights;            %- Vector of 1X4 with the weights given to the attributes of classes 1-4 based on their relative information gain.       
     P.All_objs=All_objs;              % - struct of 1X5 containing the GMM parameters of the clicks' similarity attributes
  %% import Glider metadata   
     load Glider_meta_data;   % load table with glider's IMU measurements
     P.AllData=AllData;       % glider's IMU measurements

end