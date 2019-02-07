function [Pcov_h, Pcov_m] = simulate_uplink_coverage_with_nrb(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.N_RB);
PoCoverageH = zeros(params.space_realizations , params.time_slots );
PoCoverageM = zeros(params.space_realizations , params.time_slots );
for p = 1:points
fprintf('\n')
disp(['Number of Resource Blocks: ' , num2str(params.N_RB(p))]);
disp(['Smallcell Density: ' , num2str(params.LA_B*1e6)]);
disp(['HTC  Density: ' , num2str(params.LA_H*1e6)]);
disp(['MTC Density: ' , num2str(params.LA_M)]);
disp(['MTC Activity Ratio: ' , num2str(params.rho_m)]);
disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
for m = 1:params.space_realizations;
    if(mod(m,params.space_realizations/100) == 0)
        fprintf('|');
    end
    mu_b = params.LA_B * simulation_area;
    mu_m = params.LA_M * simulation_area;
    mu_h = params.LA_H * simulation_area;
    
    N_cells = poissrnd(mu_b);
    while(N_cells < 1)
        N_cells = poissrnd(mu_b);
    end
    N_users_MTC = round(params.rho_m * poissrnd(mu_m));
    while(N_users_MTC < 1)
        N_users_MTC = round(params.rho_m * poissrnd(mu_m));
    end
    N_users_HTC = round(poissrnd(mu_h));
    while(N_users_HTC < 1)
        N_users_HTC = poissrnd(mu_h);
    end
    locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
    locations.HTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_HTC, 2);
    locations.MTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_MTC, 2);
    
    % Association
    
    [association,loads,distances] = associate_nodes(locations,params.aggregation_mode);
    N_cells_active = association.No_Active_Cells;
    RBs = randi(params.N_RB(p),1,N_users_MTC);

    association.MTC_to_BS(2,:)= RBs(1,:);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
    
    % Power Control
    server_MTC_P = params.MTC.Pmin * exp(params.SEPL.alpha .* distances.MTC_to_Server .^ params.SEPL.beta);
    server_MTC_P(server_MTC_P > params.MTC.Pmax) = 0;
    server_HTC_P = params.HTC.Pmin * exp(params.SEPL.alpha .* distances.HTC_to_Server .^ params.SEPL.beta);
    server_HTC_P(server_HTC_P > params.HTC.Pmax) = 0;
    
    % Pairing of HTC users with MTC nodes  and Generating of ordered channels
    for server = 1:N_cells_active
        HTC_MTC_Pair(server).HTCs = find( association.HTC_to_BS == server);
        HTC_MTC_Pair(server).MTCs = find( association.MTC_to_BS(1,:) == server);
        for r = 1:params.N_RB(p)
            y = [params.H(r,1) params.H(r,2)];
            ordered(r,:) = order_exp(y);
        end
    end
    
    Interfers_h = zeros(N_cells_active,1);
    for server = 1:N_cells_active
         interferer = find( association.HTC_to_BS == server,1,'first');
         if(numel(interferer) ~= 0)
             Interfers_h(server) = interferer;
         end
    end
    for t = 1:params.time_slots
        S_m = zeros(N_users_MTC,1);
        I_m = zeros(N_users_MTC,1);
        I_m_h = zeros(N_users_MTC,1);
        Isic = zeros(N_users_MTC,1);
        S_h = zeros(N_users_HTC,1);
        I_h = zeros(N_users_HTC,1);
        I_h_m = zeros(N_users_HTC,1);
        Inoma = zeros(N_users_HTC,1);
%         H_m = exprnd(1,N_cells_active,N_users_MTC);
%         H_h = exprnd(1,N_cells_active,N_users_HTC);
        H_m = params.H(m+t:N_cells_active+m+t-1 , m+t: N_users_MTC+m+t-1);
        toss = randi(10);
        H_h = params.H(m+t+toss:N_cells_active+m+t+toss-1 , m+t+toss: N_users_HTC+m+t+toss-1);

        for srvr = 1:N_cells_active
            htcs = HTC_MTC_Pair(srvr).HTCs;
            mtcs = HTC_MTC_Pair(srvr).MTCs;
            indx = mod(srvr,params.N_RB);
            H_m(srvr,mtcs) = ordered(indx+1,1);
            H_h(srvr,htcs) = ordered(indx+1,2);
        end
%         for u_m = 1:N_users_MTC
%             server_m = association.MTC_to_BS(1,u_m);
%             H_m(server_m,u_m) = ordered_h(mod(u_m,params.N_RB(p))+1,1) ;
%         end
%         for u_h = 1:N_users_HTC
%             server_h = association.HTC_to_BS(1,u_h);
%             H_h(server_h,u_h) = ordered_h(mod(u_h,params.N_RB(p))+1,2);
%         end
        
        
        R_h = H_h .* bsxfun(@times,server_HTC_P, exp(-params.SEPL.alpha .* distances.HTC_to_BS .^ params.SEPL.beta))  ;
        %R_m = H_m .* bsxfun(@times,server_MTC_P, exp(-params.SEPL.alpha .* distances.MTC_to_BS .^ params.SEPL.beta))  ;

        for user_m = 1:N_users_MTC
            server_m = association.MTC_to_BS(1,user_m);
            rb = association.MTC_to_BS(2,user_m);
            %same_ser =  (association.MTC_to_BS(1,:) == server_m);
            another_ser =  (association.MTC_to_BS(1,:) ~= server_m);
            same_rb = (association.MTC_to_BS(2,:) == rb);
            %intra_interferers = find(same_ser + same_rb == 2);
            MTC_interferers = find(another_ser + same_rb == 2);
            S_m(user_m) = params.MTC.Pmin * H_m(server_m,user_m);
            %I_intra =  sum(params.MTC.Pmin * H_m(server_m,intra_interferers));
            %distance_to_interferers = pdist2(locations.BS(server_m,:),locations.MTC(inter_interferers,:),'euclidean') ;
            distance_to_MTC_interferers = distances.MTC_to_BS(server_m,MTC_interferers);
            I_inter = sum( server_MTC_P(MTC_interferers) .* H_m(server_m,MTC_interferers).* exp(-params.SEPL.alpha .* distance_to_MTC_interferers .^ params.SEPL.beta) );
            %I_m(i) = I_intra + I_inter;
            I_m(user_m) =  I_inter;  % special case #1 no intra interference
            %HTC_interferers = setdiff(1:N_users_HTC ,
            %connected_users{server_m});k
            if(Interfers_h(server_m) ~= 0)
                I_m_h(user_m) = sum(R_h(server_m,Interfers_h(Interfers_h~=0))) - R_h(server_m,Interfers_h(server_m));
            else
                I_m_h(user_m) = sum(R_h(server_m,Interfers_h(Interfers_h~=0)));
            end
            Isic(user_m) = H_m(server_m,user_m) * params.HTC.Pmin;
        end
        for user_h = 1:N_users_HTC
            server_h = association.HTC_to_BS(user_h);
            users_m = find(association.MTC_to_BS(1,:) == server_h); % find MTC nodes associated to this server
            user_m = users_m(1); 
            gh = H_h(server_h,user_h);
            gm = H_m(server_h,user_m);
            x = order_exp([gm gh]);
            S_h(user_h) = x(2)/gh*R_h(server_h,user_h);
            if(Interfers_h(server_h) ~= 0)
                I_h(user_h) = sum(R_h(server_h,Interfers_h(Interfers_h~=0))) -  R_h(server_h,Interfers_h(server_h));
            else
                I_h(user_h) = sum(R_h(server_h,Interfers_h(Interfers_h~=0))) ;
            end 
            %I_h_(user_h) = sum(R_h(server_h,:)) -  sum(R_h(server_h,connected_users{server_h})) ;  % sum all signals received at this server from outside its coverage area.
            if(numel(users_m) > 0)
                %I_h_m(user_h) =  I_m(users_m(1))   ;     % the interfernce on this user from MTC nodes is the interference on any MTC node (say 1)
                I_h_m(user_h) =  mean(I_m(users_m) )  ; 
                Inoma(user_h) = x(1)/gm*S_m(users_m(1));
                %Inoma(user_h) = mean(S_m(users_m));
            end
        end
%         for user_m = 1:N_users_MTC
%             server_m = association.MTC_to_BS(1,user_m);
%             HTC_interferers = setdiff(1:N_users_HTC , connected_users{server_m});
%             %distances_to_HTC_interferers = pdist2(locations.BS(server_m,:),locations.HTC(HTC_interferers,:),'euclidean') ;
%             distances_to_HTC_interferers = distances.HTC_to_BS(server_m,HTC_interferers);
%             R_h = H_h(server_m,HTC_interferers) .* bsxfun(@times,server_HTC_P(HTC_interferers), exp(-params.SEPL.alpha .* distances_to_HTC_interferers.^ params.SEPL.beta))  ;
%             I_m_h(user_m) = sum(R_h);
%             Isic = H_m(server_m,user_m) * params.HTC.Pmin;
%         end
         SINR_h = S_h ./ (I_h + I_h_m + Inoma +  params.No);   % Noise here is per only one RB !!!!!!
         SINR_m = S_m ./ (I_m +  I_m_h + params.No);% +  +    % Noise here is per only one RB !!!!!!
         for u_m = 1:N_users_MTC
             server_m = association.MTC_to_BS(1,u_m);
             u_h = HTC_MTC_Pair(server_m).HTCs;
             if(SINR_h(u_h(1)) < params.Threshold.HTC_QOS)
                SINR_m(u_m) = 0;
             end   
         end
         
            PoCoverageH(m,t) = sum(SINR_h > params.Threshold.HTC) / N_users_HTC;
            PoCoverageM(m,t) = sum(SINR_m > params.Threshold.MTC) / N_users_MTC;

    end
    
end

    normfact = params.space_realizations * params.time_slots ;% * simulation_area;
    Pcov_h(p) = sum(sum(PoCoverageH))/ normfact
    Pcov_m(p) = sum(sum(PoCoverageM))/ normfact
    fprintf('\n');

display('\');
toc;
end
end