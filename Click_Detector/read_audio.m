function [Raw_audio,Fs,Edge_flag]=read_audio(DF,PF,Audio_name,buffer_size,buffer_index)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           July 2025
%   DESCRIPTION:

% read_audio Returns a specified segment from a given audio file  

%   INPUT:
%   > DF               - Directory of the data folder
%   > PF               - Directory of the program folder
%   > Audio_name       - Full_name.format of the audio file
%   > buffer_size      - User defined window size for analysis 
%   > buffer_index     - Start index for buffer selection 

%   OUTPUT:
%   > Raw_audio      - Vector of MX4 comprising 4 channel of M samples of a 30-second segment from a the audio file
%   > Fs             - Scalar representing the sampled frequency of the audiofile
%   > Edge_flag      - Scalar that is set to 1 if the buffer index exceeds the audio length, and is set to 0 otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Edge_flag=0;  
    cd(DF)
    Audio_info=audioinfo(Audio_name);
    Fs=Audio_info.SampleRate;

    if buffer_index*buffer_size<=Audio_info.Duration
        disp('Reading audio...')
        Raw_audio = audioread(Audio_name,[(buffer_index-1)*buffer_size*Fs+1,buffer_index*buffer_size*Fs]);            % read recording
    else
        Raw_audio=[];
        Edge_flag=1;
    end
    cd(PF)
    
    

end