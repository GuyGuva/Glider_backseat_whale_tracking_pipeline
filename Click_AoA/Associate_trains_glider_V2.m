function Whale_output=Associate_trains_glider_V2(P,record_time,DoA_trajectories,trajectories,id_j_ToAs,test,Fs,Plot_flag)
  
%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2025
%   DESCRIPTION:

% This function assigns click trains into sources. This assignment involves 
% grouping trains of similar distributions in terms of the slant delay while
% taking into account possible silent periods between click trains due to 
% missed detections, foraging, or breathing cycles
% The function outputs a struct containing the arrival times and slant
% delays of all clicks of each identified sorce whale

%   INPUT:
%   > trajectories       - Sturct containing the indices of click trains
%   > id_j_ToAs          - Sturct containing the clicks' arrival times of click trains
%   > DoA_trajectories   - Cell array containing the clicks' arrival angles within the formed trains
%   > test               - Vector reprsenting a band-passed signal
%   > t1                 - Scalar indicating the start time (in datetime format) of the audio 
%   > Fs                 - Scalar: smpale rate of the signal (test)
%   > P                  - Struct containing the detector parameters
%   > Plot_flag          - Flag for visualizing results

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Whale(1).ToAs=[];
    Whale(1).Elevation=[];
    Whale(1).Azimuth=[];
    Whale(1).Pks=[];
    trajectories = trajectories(~cellfun(@isempty, trajectories));
    Nd=length(trajectories);
    if Nd>0
        L_tot = 1e9*ones(Nd,Nd);
        for i=1:Nd-1
            locs_i=cell2mat(id_j_ToAs(trajectories{i})');
            t_i_last=locs_i(end);                
            ref_ToAs=cell2mat(DoA_trajectories(i));
            mu_i=mean(ref_ToAs);
            Sigma_i=eye(2).*std(ref_ToAs);
            for j=i+1:Nd  
                locs_j=cell2mat(id_j_ToAs(trajectories{j})');
                t_j_first=locs_j(1);
                ref_ToAs=cell2mat(DoA_trajectories(j));
                mu_j=median(ref_ToAs);
                Sigma_j=eye(2).*std(ref_ToAs);
    
                if t_j_first-t_i_last>P.ITI_min && t_j_first-t_i_last<P.ITI_max                          
                     L_tot(i,j)=bhattacharyya_gaussian(mu_i', Sigma_i, mu_j', Sigma_j);                                       
                end  
    
           end
        end        
        [C,~] = do_assignment (L_tot,P.lone_p, 1e9);
    
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
    
        individual_trains_idx=find(chain_nr==0);
        Max_idx=max(chain_nr);
        for ii=1:length(individual_trains_idx)
          L_train=length(cell2mat(id_j_ToAs(cell2mat(trajectories(individual_trains_idx(ii))))'));
          if L_train>P.min_numer_of_clicks_per_whale
              Max_idx=Max_idx+1;
              chain_nr(individual_trains_idx(ii))=Max_idx;
          end
        end
     
         Det={};
         c=0;
         for i=1:max(chain_nr)
             ind=find(chain_nr==i); 
             if ~isempty(ind)
                 c=c+1;
                 Det(c)={ind};
             end
         end
         
         if ~isempty(Det)         
             co={'r*','g*','b*','c*','m*','k*','y*','r+','g+','b+','c+','m+','k+','y+','ro','go','bo','co','mo','ko','yo'};        
             if Plot_flag.Click_AoA_final_time_series
                figure;
                subplot(3,1,1); hold off;
                t_test=(0:1/Fs: (1/Fs)*(length(test)-1))';
                plot(t_test,test);  hold on; grid on;
                subplot(3,1,2); hold off; 
                plot(0,0,'.'); hold on; grid on;  xlabel('Time [sec]'); ylabel('Elevation [\circ]'); 
                subplot(3,1,3); hold off; 
                plot(0,0,'.'); hold on; grid on; xlabel('Time [sec]'); ylabel('Azimuth [\circ]'); 
             end
             cc=1;
             legendInfo{cc} = '';
             Det = Det(~cellfun(@isempty, Det));
             Det_EL=[];
             for i=1:length(Det)
                whale1=Det{i}';
                sol=cell2mat(id_j_ToAs(cell2mat(trajectories(whale1)))');
                if length(sol)>P.min_numer_of_clicks_per_whale
                    cc=cc+1;
                    DoA_whale=DoA_trajectories(whale1); 
                    DoA_whale_Matrix = vertcat(DoA_whale{:});
                    Whale(cc-1).ToAs=sol;
                    Whale(cc-1).Elevation=DoA_whale_Matrix(:,1)';
                    Whale(cc-1).Azimuth=DoA_whale_Matrix(:,2)';
                    Pks=Peaks_extract(test,sol,Fs);
                    Whale(cc-1).Pks=Pks;
                    if Plot_flag.Click_AoA_final_time_series
                        subplot(3,1,2); 
                        dotH=plot(sol,DoA_whale_Matrix(:,1)',co{cc-1},'LineWidth',2);  hold on;
                        pause(0.005);
                        set(dotH, 'XData', sol, 'YData', DoA_whale_Matrix(:,1)'); hold on;       
                        legendInfo{cc} = ['Whale ' num2str(cc-1)];
                        % legend(legendInfo)
    
                        subplot(3,1,3); 
                        dotH=plot(sol,DoA_whale_Matrix(:,2)',co{cc-1},'LineWidth',2);  hold on;
                        pause(0.005);
                        set(dotH, 'XData', sol, 'YData', DoA_whale_Matrix(:,2)'); hold on;       
                        legendInfo{cc} = ['Whale ' num2str(cc-1)];
                        % legend(legendInfo)
    
                        subplot(3,1,1); hold on;
                        dotH3=plot(sol',Pks,co{cc-1},'LineWidth',2);  hold on;
                        pause(0.005);
                        set(dotH3, 'XData', sol, 'YData', Pks); hold on;
                        % legend(legendInfo)
                    end
                else
                    Det_EL=[Det_EL i];
                    Pks_mean(i)=0;
                    Leg_txt_raw(i+1)={''};
                    Leg_txt_tdoa(i+1)={''};  
                end
                % pause;        
            end
    
             % legend(legendInfo)
             if ~isempty(Det_EL)
               Det(Det_EL)={[]};
             end
             Det = Det(~cellfun(@isempty, Det));
         end
    end

    if Plot_flag.Click_AoA_final_time_series_smoothed
        if ~isempty(Det)
            figure;
            for Q=1:size(Whale,2)
                    sol=Whale(Q).ToAs;
                    El=movmedian(movmedian(Whale(Q).Elevation,length(sol)),length(sol));
                    Az=movmedian(movmedian(Whale(Q).Azimuth,length(sol)),length(sol));
                    subplot(2,1,2); plot(sol,Az,co{Q},'LineWidth',2); hold on; grid on; ylim([-180 180]);
                    xlabel('Time [sec]'); ylabel('Azimuth [\circ]');
                    subplot(2,1,1); plot(sol,El,co{Q},'LineWidth',2); hold on; grid on; ylim([-180 180]);
                    xlabel('Time [sec]'); ylabel('Elevation [\circ]');
            end
        end
    end


   Bearings=[];
   if Plot_flag.Click_AoA_final
        figure; subplot(1,3,1);
        t_test=(0:1/Fs: (1/Fs)*(length(test)-1))';
        plot(t_test,test);  hold on; grid on;
   end
   if ~isempty(Det)
       for Q=1:size(Whale,2)
                sol=Whale(Q).ToAs;
                El=Whale(Q).Elevation;
                Az=Whale(Q).Azimuth;
                El_smoothed=movmedian(movmedian(Whale(Q).Elevation,length(sol)),length(sol));
                Az_smoothed=movmedian(movmedian(Whale(Q).Azimuth,length(sol)),length(sol));

                if Plot_flag.Click_AoA_final
                    subplot(1,3,1);
                    plot(Whale(Q).ToAs,Whale(Q).Pks,co{Q},'LineWidth',2);
                    subplot(1,3,2)
                    plot(Az_smoothed,El_smoothed,co{Q},'LineWidth',2); hold on; grid on; axis([-180 180 -180 180])
                    ylabel('Elevation'); xlabel('Azimuth');
                end
                Whale(Q).Elevation=El;
                Whale(Q).Azimuth=Az;
                Whale(Q).Elevation_smoothed=El_smoothed;
                Whale(Q).Azimuth_smoothed=Az_smoothed;
                Locs_datetime=record_time+seconds(sol);
                Orientation=Extract_glider_orientation(P.AllData,Locs_datetime);
                Whale(Q).heading_at_clicks=Orientation.heading_at_clicks;
                Whale(Q).pitch_at_clicks=Orientation.pitch_at_clicks;
                Whale(Q).roll_at_clicks=Orientation.roll_at_clicks;
                Whale(Q).Lat_at_clicks=Orientation.Lat_at_clicks;
                Whale(Q).Lon_at_clicks=Orientation.Lon_at_clicks;
                Bearings=[Bearings mean(Whale(Q).Azimuth_smoothed)];
        end
       if Plot_flag.Click_AoA_final
            subplot(1,3,3);
            visualize_on_map(1e3,Bearings,co,Whale(Q).Lat_at_clicks(1),Whale(Q).Lon_at_clicks(1),mean(Whale(Q).heading_at_clicks)); hold on;
       end
       for Q=1:size(Whale,2)
              Whale_output(Q).click_arrival_times_datetime=record_time+seconds(Whale(Q).ToAs);
              Whale_output(Q).click_arrival_times_sec=Whale(Q).ToAs;
              Whale_output(Q).click_amplitudes=Whale(Q).Pks;
              Whale_output(Q).glider_latitude=Whale(Q).Lat_at_clicks;
              Whale_output(Q).glider_longitude=Whale(Q).Lon_at_clicks;
              Whale_output(Q).glider_heading=Whale(Q).heading_at_clicks;
              Whale_output(Q).whale_bearing_from_glider=Whale(Q).Azimuth;
              Whale_output(Q).whale_smoothed_bearing_from_glider=Whale(Q).Azimuth_smoothed;
              Whale_output(Q).whale_elevation_from_glider=Whale(Q).Elevation;
              Whale_output(Q).whale_smoothed_elevation_from_glider=Whale(Q).Elevation_smoothed;              
        end
    
       Power=zeros(1,length(Whale_output));
       for i=1:length(Whale_output)
            Power(i)=mean(Whale_output(i).click_amplitudes);
            Locs_datetime=Whale_output(i).click_arrival_times_datetime;
            TT_depth = table2timetable(P.AllData, 'RowTimes','Timestamp');
            TT_matched_table = retime(TT_depth, Locs_datetime, 'nearest');
            Whale_output(i).Glider_depth=TT_matched_table.Depth;
       end
    
       lg(1)={''};
       for i=1:length(Whale_output)
            Power(i)=mean(Whale_output(i).click_amplitudes);
            NOC(i)=length(Whale_output(i).click_amplitudes);
            Rank(i)=Power(i).*NOC(i);
            Whale_output(i).Rank=Rank(i);
            lg(i+1)={num2str(round(Rank(i),2))};
        end
       if Plot_flag.Click_AoA_final
          subplot(1,3,1); legend(lg)
       end
   else
       Whale_output=[];
   end

end

