function Orientation=Extract_glider_orientation(AllData,Locs_datetime)
      
        %% 1) Ensure AllData.Timestamp is datetime with full date and time
        AllData.Timestamp = datetime(AllData.Timestamp, ...
            'InputFormat','dd/MM/yyyy HH:mm:ss');
        
        %% 2) Ensure Locs_datetime is also datetime (with the same date reference)
        % Example:
         Locs_datetime = datetime(Locs_datetime, 'InputFormat','dd/MM/yyyy HH:mm:ss');
        
        %% 3) Convert AllData into a timetable indexed by Timestamp
        TT = table2timetable(AllData, 'RowTimes', 'Timestamp');
        
        %% 4) Resample at click times using nearest-neighbor (full date + time)
        TT_clicks = retime(TT, Locs_datetime, 'nearest');
        
        %% (Optional) The actual matched timestamps:
        matched_times = TT_clicks.Timestamp;
        
        % 4) Extract heading, pitch, roll at click times
        Orientation.heading_at_clicks = TT_clicks.Heading;
        Orientation.pitch_at_clicks   = TT_clicks.Pitch;
        Orientation.roll_at_clicks    = TT_clicks.Roll;
        Orientation.Lat_at_clicks=TT_clicks.Lat;
        Orientation.Lon_at_clicks=TT_clicks.Lon;
        Orientation.declination_at_clicks=TT_clicks.Declination;
        

end