function connect_nodes_to_servers(locations,association,mode)
switch(mode)
    case 'AGGREGATION'
        for i = 1:size(association.M2M_to_H2H,1);
            x = [locations.M2M(association.M2M_to_H2H(i,1),1) , locations.H2H(association.M2M_to_H2H(i,2),1) ];
            y = [locations.M2M(association.M2M_to_H2H(i,1),2) , locations.H2H(association.M2M_to_H2H(i,2),2) ];
            line(x ,y,'Color','green','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.M2M_to_CH,1);
            x = [locations.M2M(association.M2M_to_CH(i,1),1) , locations.CH(association.M2M_to_CH(i,2),1) ];
            y = [locations.M2M(association.M2M_to_CH(i,1),2) , locations.CH(association.M2M_to_CH(i,2),2) ];
            line(x ,y,'Color','magenta','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.H2H_to_BS,2);
            x = [locations.H2H(i,1) , locations.BS(association.H2H_to_BS(i),1) ];
            y = [locations.H2H(i,2) , locations.BS(association.H2H_to_BS(i),2) ];
            line(x ,y,'Color','green','LineStyle','-','LineWidth',3)
        end
        
        for i = 1:size(association.CH_to_BS,2);
            x = [locations.CH(i,1) , locations.BS(association.CH_to_BS(i),1) ];
            y = [locations.CH(i,2) , locations.BS(association.CH_to_BS(i),2) ];
            line(x ,y,'Color','magenta','LineStyle','-','LineWidth',3)
        end
    case 'C2A'
        for i = 1:size(association.MTC_to_BS,2);
            server_m = association.MTC_to_BS(i);
            x = [locations.MTC(i,1) , locations.BS(association.Active_Servers(server_m),1) ];
            y = [locations.MTC(i,2) , locations.BS(association.Active_Servers(server_m),2) ];
            line(x ,y,'Color','red','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.HTC_to_BS,2);
            server_h = association.HTC_to_BS(i);
            x = [locations.HTC(i,1) , locations.BS(association.Active_Servers(server_h),1) ];
            y = [locations.HTC(i,2) , locations.BS(association.Active_Servers(server_h),2) ];
            line(x ,y,'Color','green','LineStyle','-','LineWidth',3)
        end
end
end