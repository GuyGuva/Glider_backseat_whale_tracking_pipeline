function [trajectories,id_j_ToAs,DoA_trajectories]=run_train_DoA_version(Detections,test,Fs,P,Plot_flag)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2025
%   DESCRIPTION:

% This function combines assigned clicks from consecutive time buffers into complete click trains. 
% This is performed by tracking the click subsets via a MAP approach
% The function outputs structs containing the formed click trains and their
% clicks' arrival times.

%   INPUT:
%   > Detections         - Sturct containing the clicks' attributes of each source in each buffer
%   > test               - Vector reprsenting a band-passed signal
%   > Fs                 - Scalar: smpale rate of the signal (test)
%   > P                  - Struct containing the detector parameters
%   > Plot_flag          - Flag for visualizing results

%   OUTPUT:
%   > trajectories     - Cell array containing the indices of the formed click trains
%   > id_j_ToAs        - Cell array containing the clicks' arrival times within the formed trains
%   > DoA_trajectories - Cell array containing the clicks' arrival angles within the formed trains

% The first two outputs are structured such that
% sol=cell2mat(id_j_ToAs(trajectories{j})'); is a vector containing the
% arrival times of all clicks in train with index j
% DoA_trajectories={Elevation,Azimuth};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    l_prev=0;
    transition_arcs=[];
    detection_arcs=[];
    co_1=0;
    co_2=0;
    for i=1:length(Detections)    
        if i==1
            C_i_en=0;
        else
            Pen=length(Detections(i).ICI)/(length(Detections(i).ICI)+length(Detections(i-1).ICI));
            if Pen==0
                Pen=eps;
            end
            C_i_en=-log(Pen);
        end
        if i<length(Detections)
            Pex=length(Detections(i).ICI)/(length(Detections(i).ICI)+length(Detections(i+1).ICI));
            if Pex==0
                Pex=eps;
            end        
           C_i_ex=-log(Pex);
        else
           C_i_ex=0;
        end
    

        C_i=log(0.01/0.99);
    
        for j=1:length(Detections(i).ICI)
            id_j=l_prev+j;       
            co_1=co_1+1;
            id_j_ToAs(co_1)=Detections(i).ToAs(j);
            DoAs_tmp=Detections(i).DoA_all_clicks{j};

            id_j_Azimuth(co_1)={DoAs_tmp(:,1)'};
            id_j_Elevation(co_1)={DoAs_tmp(:,2)'};

            detection_arcs(co_1,:)=[id_j C_i_en C_i_ex C_i];
            Pkk_current=median(Detections(i).Pkk{j});
            DoA_current=Detections(i).DoA{j};
            ICI_current=median(Detections(i).ICI{j});
    
            for k=1:length(Detections(i+1).ICI)
                id_k=l_prev+length(Detections(i).ICI)+k;
                Pkk_next=median(Detections(i+1).Pkk{k});
                DoA_next=Detections(i+1).DoA{k};
                ICI_next=median(Detections(i+1).ICI{k});

                ToA_prev=sort(Detections(i).ToAs{j});
                ToA_follow=sort(Detections(i+1).ToAs{k});
                if isempty(ToA_prev) | isempty(ToA_follow)
                    ICI_constraint=true;
                else
                    ICI_constraint=ToA_follow(1)-ToA_prev(end)>2;
                end
                alpha_Pkk=log(Pkk_next/Pkk_current);

                alpha_DoA= pdist([DoA_current ; DoA_next]);

                if abs(alpha_DoA)>P.max_AoA_change
                    DOA_constraint=1;
                else
                    DOA_constraint=0;
                end

                alpha_ICI=log(ICI_next/ICI_current);
                
                Features=[alpha_ICI alpha_Pkk alpha_DoA];
                L(j,k)=-log(Buffer_pairing_likelihood_glider(Features,P.Buffer_Params));
                if ICI_constraint || DOA_constraint
                    L(j,k)=-log(eps);
                end
    
                co_2=co_2+1;
                transition_arcs(co_2,:)=[id_j id_k L(j,k)];
                          
            end        
        end
        l_prev=l_prev+length(Detections(i).ICI);
    
    end
    
    if isempty(transition_arcs)
        c=0;
        for x=1:length(Detections)
            tmp=Detections(x).ToAs;
            if ~isempty(tmp)
            c=c+1;
            trajectories(c)={c};
            end
        end
        trajectories = trajectories';
    else    

        writematrix(detection_arcs,'detection_arcs.txt')
        writematrix(transition_arcs,'transition_arcs.txt')
    
    
        %% load affinity scores
        d_m = dlmread('detection_arcs.txt');
        t_m = dlmread('transition_arcs.txt');
        
        
        %% call mcc4mot

        [trajectories, ~] = mcc4mot(d_m,t_m);   

        %% Remove empty cells
        Remove_idx=[];
        for j=1:length(trajectories)
              if isempty(cell2mat(id_j_ToAs(trajectories{j})')')
                  Remove_idx=[Remove_idx j];
              end
        end
        trajectories(Remove_idx)=[];
        
        %% Results visualization:
   
        if Plot_flag.Click_AoA_trains_formation
            figure;
            subplot(3,1,1); hold off;
            t_test=(0:1/Fs: (1/Fs)*(length(test)-1))';
            plot(t_test,test);  hold on; grid on;
            subplot(3,1,2); hold off; 
            plot(0,0,'.'); hold on; grid on;  xlabel('Time [sec]'); ylabel('Elevation [\circ]');
            subplot(3,1,3); hold off; 
            plot(0,0,'.'); hold on; grid on; xlabel('Time [sec]'); ylabel('Azimuth [\circ]');
            subplot(3,1,2); 
            xlabel('Time [sec]'); ylabel('Elevation [\circ]'); 
        end
        for j=1:length(trajectories)
            sol=cell2mat(id_j_ToAs(trajectories{j})');
            ToAs_groups=id_j_ToAs(trajectories{j});
            Elevation=cell2mat(id_j_Elevation(trajectories{j}));
            Azimuth=cell2mat(id_j_Azimuth(trajectories{j}));

            Trains(j).Elevation=Elevation;
            Trains(j).Azimuth=Azimuth;
            Trains(j).DoA=[Elevation' Azimuth'];

            if Plot_flag.Click_AoA_trains_formation
                subplot(3,1,2);
                dotH=plot(sol,Elevation,'*','LineWidth',2);  hold on;
                pause(0.01);
                set(dotH, 'XData', sol, 'YData', Elevation); hold on;   
                subplot(3,1,3);
                dotH=plot(sol,Azimuth,'*','LineWidth',2);  hold on;
                pause(0.01);
                set(dotH, 'XData', sol, 'YData', Azimuth); hold on;    
                subplot(3,1,1); hold on;
                Pks=Peaks_extract(test,sol,Fs);  
                Pks_mean(j)=1e3*mean(Pks);
                dotH3=plot(sol',Pks,'*','LineWidth',2);  hold on;
                pause(0.01);
                set(dotH3, 'XData', sol, 'YData', Pks); hold on;
            end
            clear Elevation;
            clear Azimuth;           
        end
        init=[];
        for j=1:length(trajectories)
            sol=cell2mat(id_j_ToAs(trajectories{j})');
            init(j)=sol(1);
        end
        [~,I]=sort(init);
        trajectories=trajectories(I);
        for i=1:length(I)
            DoA_trajectories(i)={Trains(I(i)).DoA};
        end
    end

end