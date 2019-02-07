function [association , loads  , distances] = associate_nodes(locations,mode)
switch(mode)
    case 'AGGREGATION'
        distances_HTC_to_BS = pdist2(locations.BS,locations.HTC,'euclidean') ;
        [~ , server_HTC] = min(distances_HTC_to_BS);
        association.HTC_to_BS = server_HTC;
        distances.HTC_to_BS = distances_HTC_to_BS;
        
        distances_CH_to_BS = pdist2(locations.BS,locations.CH,'euclidean') ;
        [~ , server_ch] = min(distances_CH_to_BS);
        association.CH_to_BS = server_ch;
        distances.HTC_to_BS = distances_CH_to_BS;
        
        aggregators = [locations.HTC ; locations.CH];
        N_HTC = size(locations.HTC ,1);
        N_CH = size(locations.CH ,1);
        distances_MTC_to_HTC = pdist2(aggregators,locations.MTC,'euclidean') ;
        [~ , server_MTC] = min(distances_MTC_to_HTC);
        
        IMTC_HTC = find(server_MTC <= N_HTC);
        IMTC_ch = find(server_MTC > N_HTC);
        IHTC = server_MTC(IMTC_HTC);
        Ich = server_MTC(IMTC_ch) - N_HTC  ;
        
        association.MTC_to_HTC = [IMTC_HTC' IHTC'];
        association.MTC_to_CH = [IMTC_ch' Ich'];
        
        for i = 1:N_HTC
            loads.HTC(i) = sum(association.MTC_to_HTC(:,2) == i);
        end
        
        for i = 1:N_CH
            loads.CH(i) = sum(association.MTC_to_CH(:,2) == i);
        end
        
        for i = 1:size(locations.BS ,1)
            loads.UE_HTC(i) = sum(association.HTC_to_BS == i);
            loads.UE_CH(i) = sum(association.CH_to_BS == i);
        end
        
    case 'C2A'    
        distances_HTC_to_BS = pdist2(locations.BS,locations.HTC,'euclidean') ;
        [distance_to_server_HTC , server_HTC] = min(distances_HTC_to_BS);
        
        [active_servers, ~, new_servers] = unique(server_HTC);
        distances_HTC_to_active_BS = distances_HTC_to_BS(active_servers,:);
        distances.HTC_to_BS = distances_HTC_to_active_BS;
        distances.HTC_to_Server = distance_to_server_HTC;
        association.HTC_to_BS = new_servers';
        association.Active_Servers = active_servers;
        association.No_Active_Cells = numel(active_servers);
       
        
        distances_MTC_to_BS = pdist2(locations.BS(active_servers,:),locations.MTC,'euclidean') ;
        [distance_to_server_MTC , server_MTC] = min(distances_MTC_to_BS);
        association.MTC_to_BS = server_MTC;
        distances.MTC_to_BS = distances_MTC_to_BS;
        distances.MTC_to_Server = distance_to_server_MTC;
        
%         N_BS = size(locations.BS ,1);
%         for i = 1:N_BS
%             loads.BS(i,1) = sum(association.MTC_to_BS == i);
%             loads.BS(i,2) = sum(association.HTC_to_BS == i);
%         end
        
    case 'C2C'
        
        distances_HTC_to_BS = pdist2(locations.BS,locations.HTC,'euclidean') ;
        [~ , server_HTC] = min(distances_HTC_to_BS);
        distances_MTC_to_BS = pdist2(locations.BS,locations.MTC,'euclidean') ;
        [~ , server_MTC] = min(distances_MTC_to_BS);

        active_servers_all = [server_HTC server_MTC];
        [active_servers, ~, ~] = unique(active_servers_all);
        [active_servers_HTC, ~, ~] = unique(server_HTC);
        [active_servers_MTC, ~, ~] = unique(server_MTC);
        association.No_Active_Cells = numel(active_servers);
        association.No_Active_Cells_HTC = numel(active_servers_HTC);
        association.No_Active_Cells_MTC = numel(active_servers_MTC);
        
        distances_HTC_to_BS = pdist2(locations.BS(active_servers,:),locations.HTC,'euclidean') ;
        [distance_to_server_HTC , server_HTC] = min(distances_HTC_to_BS);
        association.HTC_to_BS = server_HTC;
        distances.HTC_to_BS = distances_HTC_to_BS;
        distances.HTC_to_Server = distance_to_server_HTC;
        
        distances_MTC_to_BS = pdist2(locations.BS(active_servers,:),locations.MTC,'euclidean') ;
        [distance_to_server_MTC , server_MTC] = min(distances_MTC_to_BS); 
        association.MTC_to_BS = server_MTC;
        distances.MTC_to_BS = distances_MTC_to_BS;
        distances.MTC_to_Server = distance_to_server_MTC;
%         N_BS = size(locations.BS ,1);
%         for i = 1:N_BS
%             loads.BS(i,1) = sum(association.MTC_to_BS == i);
%             loads.BS(i,2) = sum(association.HTC_to_BS == i);
%         end
end
loads = 0;
end