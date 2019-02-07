function [R_oma R_noma] = gen_initial_results(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_M);
%PoCoverage = zeros(points,  params.space_realizations , params.time_slots );
Sum_Rate_OMA = zeros(params.space_realizations , params.time_slots );
Sum_Rate_NOMA = zeros(params.space_realizations , params.time_slots );
normfact = params.space_realizations * params.time_slots * simulation_area * 1e-6;

for p = 1:points
    fprintf('\n')
    disp(['Smallcell Density: ' , num2str(params.LA_B)]);
    disp(['MTC Density: ' , num2str(params.LA_M(p))]);
    disp(['HTC Density: ' , num2str(params.LA_H)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp(['Power Control Parameters (HTC): ' ,'Pmax=', num2str(params.HTC.Pmax) ,' Pmin=',num2str(params.HTC.Pmin_dBm)]);
    disp(['Power Control Parameters (MTC): ' ,'Pmax=', num2str(params.MTC.Pmax) ,' Pmin=',num2str(params.MTC.Pmin_dBm)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        mu_h = params.LA_H * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M(p) * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_MTC = round(params.rho_m * poissrnd(mu_m));
        N_users_HTC = round(poissrnd(mu_h));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.HTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_HTC, 2);
        locations.MTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_MTC, 2);
        
        % Association
        distances_HTC_to_BS = pdist2(locations.BS,locations.HTC,'euclidean') ;
        [distance_to_server_HTC , server_HTC] = min(distances_HTC_to_BS);
        ActiveBS = unique(server_HTC);
        N_active_cells = numel(ActiveBS);
        
        locations.ActiveBS = locations.BS(ActiveBS,:);
        distances_HTC_to_Active_BS = pdist2(locations.ActiveBS,locations.HTC,'euclidean') ;
        distances_MTC_to_Active_BS = pdist2(locations.ActiveBS,locations.MTC,'euclidean') ;
        [distance_to_server_MTC , servers_mtc] = min(distances_MTC_to_Active_BS);
        
        server_MTC = ActiveBS(servers_mtc);
        for i = 1:N_cells
            loads(i,1) = sum(server_HTC == i);
            loads(i,2) = sum(server_MTC == i);
        end
        % OMA resource allocation (base line)
        Resource_Allocation_HTC_OMA = zeros(N_cells,params.N_RB);
        Resource_Allocation_HTC_OMA_Pointer = ones(N_cells,1);
        for HTC_user = 1:N_users_HTC
            HTC_user_server = server_HTC(HTC_user);
            pointer = Resource_Allocation_HTC_OMA_Pointer(HTC_user_server);
            allocation_range = pointer:pointer+params.N_RBs_HTC-1;
            Resource_Allocation_HTC_OMA(HTC_user_server,allocation_range) = HTC_user*ones(1,params.N_RBs_HTC) ;
            Resource_Allocation_HTC_OMA_Pointer(HTC_user_server) = Resource_Allocation_HTC_OMA_Pointer(HTC_user_server) + params.N_RBs_HTC;
        end
        
        Resource_Allocation_MTC_OMA = zeros(N_cells,params.N_RB);
        Resource_Allocation_MTC_OMA_Pointer = Resource_Allocation_HTC_OMA_Pointer;
        for MTC_user = 1:N_users_MTC
            MTC_user_server = server_MTC(MTC_user);
            pointer = Resource_Allocation_MTC_OMA_Pointer(MTC_user_server);
            allocation_range = pointer:pointer+params.N_RBs_MTC-1;
            Resource_Allocation_MTC_OMA(MTC_user_server,allocation_range) = MTC_user*ones(1,params.N_RBs_MTC) ;
            Resource_Allocation_MTC_OMA_Pointer(MTC_user_server) = Resource_Allocation_MTC_OMA_Pointer(MTC_user_server) + params.N_RBs_MTC;
        end
        
        % NOMA resource allocation
        Resource_Allocation_HTC_NOMA = zeros(N_cells,params.N_RB);
        Resource_Allocation_HTC_NOMA_Pointer = ones(N_cells,1);
        for HTC_user = 1:N_users_HTC
            HTC_user_server = server_HTC(HTC_user);
            pointer = Resource_Allocation_HTC_NOMA_Pointer(HTC_user_server);
            allocation_range = pointer:pointer+params.N_RBs_HTC-1;
            Resource_Allocation_HTC_NOMA(HTC_user_server,allocation_range) = HTC_user*ones(1,params.N_RBs_HTC) ;
            Resource_Allocation_HTC_NOMA_Pointer(HTC_user_server) = Resource_Allocation_HTC_NOMA_Pointer(HTC_user_server) + params.N_RBs_HTC;
        end
        
        Resource_Allocation_MTC_NOMA = zeros(N_cells,params.N_RB);
        Resource_Allocation_MTC_NOMA_Pointer = ones(N_cells,1);
        
        for MTC_user = 1:N_users_MTC
            MTC_user_server = server_MTC(MTC_user);
            pointer = Resource_Allocation_MTC_NOMA_Pointer(MTC_user_server);
            allocation_range = pointer : pointer + params.N_RBs_MTC-1;
            Resource_Allocation_MTC_NOMA(MTC_user_server,allocation_range) = MTC_user*ones(1,params.N_RBs_MTC) ;
            Resource_Allocation_MTC_NOMA_Pointer(MTC_user_server) = Resource_Allocation_MTC_NOMA_Pointer(MTC_user_server) + params.N_RBs_MTC;
        end
        
        for t = 1:params.time_slots
            Hm = exprnd(1,N_active_cells,N_users_MTC);
            Hh = exprnd(1,N_active_cells,N_users_HTC);
            HTC_uplink_P = params.HTC.Pmin * exp(params.SEPL.alpha .* distance_to_server_HTC .^ params.SEPL.beta);
            MTC_uplink_P = params.MTC.Pmin * exp(params.SEPL.alpha .* distance_to_server_MTC .^ params.SEPL.beta);
            HTC_uplink_P(HTC_uplink_P > params.HTC.Pmax) = 0;
            MTC_uplink_P(MTC_uplink_P > params.MTC.Pmax) = 0;
            
            HTC_uplink_interfernce_signal = Hh.* bsxfun(@times,HTC_uplink_P,exp(-params.SEPL.alpha .* distances_HTC_to_Active_BS .^ params.SEPL.beta)) ;
            MTC_uplink_interfernce_signal = Hm.* bsxfun(@times,MTC_uplink_P,exp(-params.SEPL.alpha .* distances_MTC_to_Active_BS .^ params.SEPL.beta)) ;
            % OMA rate calculations
            HTC_Interference_Map_OMA = zeros(N_cells,params.N_RB);
            MTC_Interference_Map_OMA = zeros(N_cells,params.N_RB);
            for r = 1:params.N_RB
                HTC_tagged = nonzeros(Resource_Allocation_HTC_OMA(:,r));
                MTC_tagged = nonzeros(Resource_Allocation_MTC_OMA(:,r));
                for c = 1:N_active_cells
                    cell = ActiveBS(c);
                    HTC_user = Resource_Allocation_HTC_OMA(cell,r);
                    MTC_user = Resource_Allocation_MTC_OMA(cell,r);
                    total_received_signal = sum(HTC_uplink_interfernce_signal(c,HTC_tagged)) + sum(MTC_uplink_interfernce_signal(c,MTC_tagged));
                    if(HTC_user ~= 0)
                        HTC_Interference_Map_OMA(cell,r) = total_received_signal - params.HTC.Pmin*Hh(c,HTC_user);
                    end
                    if(MTC_user ~= 0)
                        MTC_Interference_Map_OMA(cell,r) = total_received_signal - params.MTC.Pmin*Hm(c,MTC_user);
                    end
                end
            end
            
            HTC_Interference_Map_OMA(HTC_Interference_Map_OMA == 0) = inf;
            MTC_Interference_Map_OMA(MTC_Interference_Map_OMA == 0) = inf;
            SINR_HTC_OMA = params.HTC.Pmin ./ (HTC_Interference_Map_OMA + params.No) ;
            SINR_MTC_OMA = params.MTC.Pmin ./ (MTC_Interference_Map_OMA + params.No) ;
            
            % NOMA rate calculations
            HTC_Interference_Map_NOMA = zeros(N_cells,params.N_RB);
            MTC_Interference_Map_NOMA = zeros(N_cells,params.N_RB);
            for r = 1:params.N_RB
                HTC_tagged = nonzeros(Resource_Allocation_HTC_NOMA(:,r));
                MTC_tagged = nonzeros(Resource_Allocation_MTC_NOMA(:,r));
                for c = 1:N_active_cells
                    cell = ActiveBS(c);
                    HTC_user = Resource_Allocation_HTC_NOMA(cell,r);
                    MTC_user = Resource_Allocation_MTC_NOMA(cell,r);
                    total_received_signal = sum(HTC_uplink_interfernce_signal(c,HTC_tagged)) + sum(MTC_uplink_interfernce_signal(c,MTC_tagged));
                    if(HTC_user ~= 0)
                        HTC_Interference_Map_NOMA(cell,r) = total_received_signal - params.HTC.Pmin*Hh(c,HTC_user);
                    end
                    if(MTC_user ~= 0)
                        MTC_Interference_Map_NOMA(cell,r) = total_received_signal - params.MTC.Pmin*Hm(c,MTC_user);
                    end
                end
            end
            
            HTC_Interference_Map_NOMA(HTC_Interference_Map_NOMA == 0) = inf;
            MTC_Interference_Map_NOMA(MTC_Interference_Map_NOMA == 0) = inf;
            SINR_HTC_NOMA = params.HTC.Pmin ./ (HTC_Interference_Map_NOMA + params.No) ;
            SINR_MTC_NOMA = params.MTC.Pmin ./ (MTC_Interference_Map_NOMA + params.No) ;
            
            Rate_HTC_OMA = params.RB_BW * log(1 + SINR_HTC_OMA);
            Rate_MTC_OMA = params.RB_BW * log(1 + SINR_MTC_OMA);
            Rate_HTC_NOMA = params.RB_BW * log(1 + SINR_HTC_NOMA);
            Rate_MTC_NOMA = params.RB_BW * log(1 + SINR_MTC_NOMA);
            
            Sum_Rate_OMA(m,t) = sum(sum(Rate_HTC_OMA  + Rate_MTC_OMA));
            Sum_Rate_NOMA(m,t)= sum(sum(Rate_HTC_NOMA + Rate_MTC_NOMA));
        end
    end
    R_oma(p) = sum(sum(Sum_Rate_OMA)) / normfact
    R_noma(p) = sum(sum(Sum_Rate_NOMA)) / normfact
end


end