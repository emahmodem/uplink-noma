function [Pcov_h, Pcov_m] = simulate_downlink_coverage_with_htc_density(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_H);
PoCoverageH = zeros(params.space_realizations , params.time_slots );
PoCoverageM = zeros(params.space_realizations , params.time_slots );

for p = 1:points
    fprintf('\n')
    %disp(['coverage Threshold: ' , num2str(params.Threshold(p))]);
    disp(['Smallcell Density: ' , num2str(params.LA_B*1e6)]);
    disp(['HTC Density: ' , num2str(params.LA_H(p)*1e6)]);
    disp(['M2M Density: ' , num2str(params.LA_M*1e6)]);
    disp(['M2M Activity Ratio: ' , num2str(params.rho_m)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        mu_m = params.LA_M * simulation_area;
        mu_h = params.LA_H(p) * simulation_area;
        
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
        
        %PA(m) = N_cells_active / N_cells;
        S_m = zeros(N_users_MTC,1);
        S_h = zeros(N_users_HTC,1);
        I_m = zeros(N_users_MTC,1);
        I_h = zeros(N_users_HTC,1);
        
        for t = 1:params.time_slots
            %             H_m = exprnd(1,N_cells_active,N_users_MTC);
            %             H_h = exprnd(1,N_cells_active,N_users_HTC);
            %                     H_m = params.H(t:N_cells_active+t-1 , t: N_users_MTC+t-1);
            %                     H_h = params.H(t:N_cells_active+t-1 , t: N_users_HTC+t-1);
            H_m = params.H(m+t:N_cells_active+m+t-1 , m+t: N_users_MTC+m+t-1);
            toss = randi(10);
            H_h = params.H(m+t+toss:N_cells_active+m+t+toss-1 , m+t+toss: N_users_HTC+m+t+toss-1);
            R_h =  params.PD *  H_h.* exp(- params.SEPL.alpha.*distances.HTC_to_BS.^params.SEPL.beta);
            R_m =  params.PD *  H_m.* exp(- params.SEPL.alpha.*distances.MTC_to_BS.^params.SEPL.beta);
            for user_h = 1:N_users_HTC
                server_h = association.HTC_to_BS(user_h);
                S_h(user_h) = params.NOMA.epsilon * R_h(server_h,user_h);
                I_h(user_h) = sum(R_h(:,user_h)) -  R_h(server_h,user_h)  + params.NOMA.theta_D * (1 - params.NOMA.epsilon) *  R_h(server_h,user_h);
            end
            SINR_h = S_h ./ (I_h + params.No);   % Noise here is per only one RB !!!!!!
            
            for user_m = 1:N_users_MTC
                server_m = association.MTC_to_BS(user_m);
                S_m(user_m) = (1 - params.NOMA.epsilon) * R_m(server_m,user_m);
                I_m(user_m) = sum(R_m(:,user_m)) - R_m(server_m,user_m) +   params.NOMA.epsilon *  R_m(server_m,user_m);
            end
            SINR_m = S_m ./ (I_m + params.No);   % Noise here is per only one RB !!!!!!
            
            PoCoverageH(m,t) = sum(SINR_h > params.Threshold.HTC) / N_users_HTC;
            PoCoverageM(m,t) = sum(SINR_m > params.Threshold.MTC) / N_users_MTC;
            
        end
        
    end
    
    normfact = params.space_realizations * params.time_slots ;% * simulation_area;
    Pcov_h(p) = sum(sum(PoCoverageH))/ normfact
    Pcov_m(p) = sum(sum(PoCoverageM))/ normfact
    fprintf('\n');
    
    
    
    toc;
    %display(strcat(num2str(round(p/points*100)) , ' %  done!' ))
    %Pa(p) = sum(PA) /  params.space_realizations
end
end