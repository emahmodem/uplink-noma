function [ output_args ] = gen_coverage_results_mMTC_NOMA_UDN(sim)
%GEN_COVERAGE_RESULTS_NOMA  Summary of this function goes here
%   Author : Mahmoud Kamel
%   Date :  10-Jun-2018 15:29:18
%   This function generates results of simulating the effect of NOMA on MTC
%   HTC traffic
%/////////////////////////////////////////////////////////////////////////
%clc; clear;  %close all;
%hold on;


% Simulation Parameters
params.simulation_area_side = [-500 500];    % simulation area side to side
params.space_realizations = 100;
params.time_slots = 10;
params.N_RB = 100;
params.Pdl = 0.5;     % in watts
params.PD = 0.5 /  params.N_RB ;  % DL transmit power per radio block
params.BW = 20e6;

params.RB_BW = 180e3; %Hz
params.No = 10^(-17.4) * params.RB_BW * 1e-3 ;
% params.N_RBs_HTC = 25;                       % each HTC user is allocated this amount of RBs in the uplink
% params.N_RBs_MTC = 1 ;                       % each MTC user is allocated this amount of RBs in the uplink

% params.SEPL.alpha = 3e-5;
% params.SEPL.beta = 2;

% params.SEPL.alpha = 0.94;
% params.SEPL.beta = 1/2;

% params.SEPL.alpha = 0.3;
% params.SEPL.beta = 2/3;

% params.SEPL.alpha = 3e-2;
% params.SEPL.beta = 1;




% NOMA Parameters
params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
params.NOMA.theta_U =0.001;

% Power Control Parameters
params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;



% Network Parameters
params.aggregation_mode = 'C2A';
params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
params.LA_M = 0.1         ;                  % MTC nodes density  (nodes/m^2)
params.rho_m = 0.3;
params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)


params.H = exprnd(1,1.1e4,4.5e4);

% Plotting Parameters
params.sim_style = {'bo','g*' ,'rx'};
params.ana_style = {'b:','g--','r-'};
params.ana_only_style = {'b:o','g*--','rx-'};
%/////////////////////////////////////////////////////////////////////////
% Simulation Switches

%sim = 'Initial_Results';
%sim = 'Visualization';


%sim = 'DL_Coverage_CoverageThreshold_C2A';
%sim = 'DL_Coverage_CoverageThreshold_C2C';

sim = 'DL_Coverage_SmallCellsDensity_C2A';
%sim = 'DL_Coverage_SmallCellsDensity_C2C';

%sim = 'DL_Coverage_HTC_Density_C2A';
%sim = 'DL_Coverage_HTC_Density_C2C';

%sim = 'DL_Coverage_MTC_Density_C2A';
%sim = 'DL_Coverage_MTC_Density_C2C';

%sim = 'DL_Coverage_NOMA_PC_C2A';
%sim = 'DL_Coverage_NOMA_PC_C2C';


%sim = 'UL_Coverage_CoverageThreshold';
switch(sim)
    case 'DL_Coverage_NOMA_PC_C2A'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_NOMA_PC_C2A');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2A';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_H = 500*1e-6 ;
        params.LA_M = 0.1 ;
        params.rho_m = 0.1;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.NOMA.epsilon = 0.05:0.05:0.45;                   % NOMA power control parmeter   \in [0,0.5]
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_noma_pc(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_noma_pc(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.NOMA.epsilon  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.NOMA.epsilon  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end 
        
        h = get(gca,'Children') ;
        deco.xlabel = 'NOMA power control ratio $\epsilon$';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.NOMA.epsilon , Pcov_m_analy(plt,:), params.ana_style{plt}, params.NOMA.epsilon , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_NOMA_PC_C2A.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_epsilon_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2a','-depsc','-r0');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
         case 'DL_Coverage_NOMA_PC_C2C'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_NOMA_PC_C2C');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_H = 500*1e-6 ;
        params.LA_M = 0.1 ;
        params.rho_m = 0.1;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.NOMA.epsilon = 0.05:0.05:0.45;                   % NOMA power control parmeter   \in [0,0.5]
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_noma_pc(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_noma_pc(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.NOMA.epsilon  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.NOMA.epsilon  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end 
        
        h = get(gca,'Children') ;
        deco.xlabel = 'NOMA power control ratio $\epsilon$';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.NOMA.epsilon , Pcov_m_analy(plt,:), params.ana_style{plt}, params.NOMA.epsilon , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_NOMA_PC_C2C.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_epsilon_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2c','-depsc','-r0');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    case 'DL_Coverage_MTC_Density_C2A'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_MTC_Density_C2A');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2A';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_H = 500*1e-6 ;
        params.rho_m = 0.1;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_mtc_density(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_mtc_density(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_M*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_M*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of MTC nodes( $\lambda_m$ nodes/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_M*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_M*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_MTC_Density_C2A.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_lam_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_lam_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_lam_c2a','-depsc','-r0');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'DL_Coverage_MTC_Density_C2C'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_MTC_Density_C2C');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_H = 500*1e-6 ;
        params.rho_m = 0.1;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_mtc_density(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_mtc_density(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_M*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_M*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of MTC nodes( $\lambda_m$ nodes/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_M*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_M*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_MTC_Density_C2C.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_lam_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_lam_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_lam_c2c','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'DL_Coverage_HTC_Density_C2A'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_HTC_Density_C2A');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2A';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.rho_m = 0.3;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.LA_H  = (50:50:500)*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_htc_density(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_htc_density(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_H*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_H*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of HTC users( $\lambda_h$ users/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_H*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_H*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_HTC_Density_C2A.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_lah_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_lah_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_lah_c2a','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'DL_Coverage_HTC_Density_C2C'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_HTC_Density_C2C');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.rho_m = 0.3;
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.LA_H  = (50:50:500)*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_htc_density(params);
            [Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_htc_density(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_H*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_H*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of HTC users( $\lambda_h$ users/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.LA_H*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_H*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_HTC_Density_C2C.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_lah_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_lah_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_lah_c2c','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'DL_Coverage_SmallCellsDensity_C2A'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_SmallCellsDensity_C2A');
        params.simulation_area_side = [-100 100];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2A';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.rho_m = 0.3;
        params.LA_M = 0.1 ;
        %params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)
        
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        %params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_H = [100 300 500]*1e-6 ;
        for i = 1:numel(LA_H)
            params.LA_H = LA_H(i);
            [~,Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_smallcells_density(params);
            [~,Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_smallcells_density(params);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_H)
            plot(params.LA_B*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_B*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_H) ;
        ind = 0;
        for i = 1:2:2*numel(LA_H)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_H(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_H(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_H)
            plot(params.LA_B*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_B*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_SmallCellsDensity_C2A.mat  params LA_H Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_las_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_las_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_las_c2a','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'DL_Coverage_SmallCellsDensity_C2C'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_SmallCellsDensity_C2C');
        %params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.rho_m = 0.3;
        params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = 0 ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = 0 ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_M = [0.01 0.05 0.1] ;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [~,Pcov_h_analy(i,:), Pcov_m_analy(i,:)] = compute_downlink_coverage_with_smallcells_density(params);
            [~,Pcov_h_simul(i,:) ,Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_smallcells_density(params);
        end
        
        % Save results
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.LA_B*1e6  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.LA_B*1e6  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.LA_B*1e6 , Pcov_m_analy(plt,:), params.ana_style{plt}, params.LA_B*1e6 , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_SmallCellsDensity_C2C.mat  params LA_M Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_las_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_las_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_las_c2c','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'DL_Coverage_CoverageThreshold_C2A'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_CoverageThreshold_C2A');
        %params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        % Network Parameters
        params.aggregation_mode = 'C2A';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_M = 0.1         ;                  % MTC nodes density  (nodes/m^2)
        params.rho_m = 0.3;
        params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:) Pcov_m_analy(i,:)] = compute_downlink_coverage_with_coverage_threshold(params);
            [Pcov_h_simul(i,:) Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_coverage_threshold(params);
        end
        
        % Save results
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.Threshold.HTC_dB  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.Threshold.HTC_dB  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.Threshold.MTC_dB , Pcov_m_analy(plt,:), params.ana_style{plt}, params.Threshold.MTC_dB , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_CoverageThreshold_C2A.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_tau_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_tau_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_tau_c2a','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'DL_Coverage_CoverageThreshold_C2C'
        tic;
        % Simulation Parameters
        disp('DL_Coverage_CoverageThreshold_C2C');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_M = 0.1         ;                  % MTC nodes density  (nodes/m^2)
        params.rho_m = 0.3;
        params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);
        
        
        
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:) Pcov_m_analy(i,:)] = compute_downlink_coverage_with_coverage_threshold(params);
            [Pcov_h_simul(i,:) Pcov_m_simul(i,:)] = simulate_downlink_coverage_with_coverage_threshold(params);
        end
        
        % Save results
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.Threshold.HTC_dB  , Pcov_h_analy(plt,:), params.ana_style{plt} ,params.Threshold.HTC_dB  ,Pcov_h_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Downlink coverage probability';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.Threshold.MTC_dB , Pcov_m_analy(plt,:), params.ana_style{plt}, params.Threshold.MTC_dB , Pcov_m_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        save DL_Coverage_CoverageThreshold_C2C.mat  params LA_B Pcov_h_analy Pcov_h_simul Pcov_m_analy Pcov_m_simul
        savefig('Downlink_Initial_Results/dl_coverage_tau_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_tau_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_tau_c2c','-depsc','-r0');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'UL_Coverage_CoverageThreshold'
        tic;
        disp('UL_Coverage_CoverageThreshold');
        %params.simulation_area_side = [-250 250];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        
        % Propagation Parameters
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 1/2;
        
        % NOMA Parameters
        params.NOMA.epsilon = 0.25;                   % NOMA power control parmeter   \in [0,0.5]
        params.NOMA.theta_D = 0.05;                        % NOMA Interfernce fraction due to SIC error propagation
        params.NOMA.theta_U =0.001;
        
        % Power Control Parameters
        params.MTC.Pmax = 2e-3 ;                    % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.MTC.Pmin_dBm = -70;                   % Min. received power at the base station in dBms > BS sensitivity
        params.MTC.Pmin = 10^(params.MTC.Pmin_dBm/10)* 1e-3;
        params.HTC.Pmax = 20e-3 ;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
        params.HTC.Pmin_dBm = -60;                   % Min. received power at the base station in watts > BS sensitivity
        params.HTC.Pmin = 10^(params.HTC.Pmin_dBm/10)* 1e-3 ;
        
        
        
        % Network Parameters
        params.aggregation_mode = 'C2C';
        params.LA_B = 1000*1e-6 ;                   % Small Cell Denisty (cells/m^2)
        params.LA_M = 0.1         ;                  % MTC nodes density  (nodes/m^2)
        params.rho_m = 0.3;
        params.LA_H = 500*1e-6 ;                  % HTC users density  (users/m^2)
        
        % Plot Independent Parameters
        params.Threshold.HTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.HTC = 10.^(params.Threshold.HTC_dB/10);
        params.Threshold.MTC_dB = (-30:0.5:10) ; % dB
        params.Threshold.MTC = 10.^(params.Threshold.MTC_dB/10);

        %         params.rho_m = 0.1;
        %         Pmin_dBm_HTC = -40 ;
        %         params.HTC.Pmin = 10.^(Pmin_dBm_HTC/10) * 1e-3;
        %         Pmin_dBm_MTC = -70 ;
        %         params.MTC.Pmin = 10.^(Pmin_dBm_MTC/10) * 1e-3;
        LA_B = [1000e-6 10000e-6 50000e-6 100000e-6];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IM = zeros(1,numel(LA_B));
        IMH = zeros(1,numel(LA_B));
        ISIC = zeros(1,numel(LA_B));
        Bh = zeros(1,numel(LA_B));
        Eph = zeros(1,numel(LA_B));
        Epm = zeros(1,numel(LA_B));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            [Pcov_h_analy(i,:) , Pcov_m_analy(i,:) ,  params, Bha(i) , Epha(i), Epma(i)] = compute_uplink_coverage_with_coverage_threshold(params);
            [Pcov_h_simul(i,:) , Pcov_m_simul(i,:) , IM(i), IMH(i),ISIC(i) , Bhs(i) , Ephs(i), Epms(i)] = simulate_uplink_coverage_with_coverage_threshold(params);
            
        end
        figure;
        subplot(1,2,1);
        %h = plot(Threshold_dB_HTC  ,Pcov_h_simul, 'r.',Threshold_dB_HTC  ,Pcov_h_analy, 'k-');
        params.sim_style = {'ko' ,'bx','gs','r*','m^'};
        params.ana_style = {'k-' ,'b-.','g--','r:','m-'};
        hold on;
        
        for plt = 1:numel(LA_B)
            plot(params.Threshold.HTC_dB  ,Pcov_h_simul(plt,:), params.sim_style{plt},params.Threshold.HTC_dB  ,Pcov_h_analy(plt,:), params.ana_style{plt});
        end
        
        h = get(gca,'Children') ;
        %h = h([4 3 2 1 ]);
        deco.xlabel = 'Coverage threshold ($\tau$)  ';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        %         deco.legend.labels = {strcat('Simulation $(\lambda_s =', num2str(LA_B(1)*1e6),'$)') ,...
        %                               strcat('Analysis $(\lambda_s =', num2str(LA_B(1)*1e6),'$)'),...
        %                               strcat('Simulation $(\lambda_s =', num2str(LA_B(2)*1e6),'$)'),...
        %                               strcat('Analysis $(\lambda_s =', num2str(LA_B(2)*1e6),'$)'),...
        %                             };
        deco.legend.location = 'northeast';
        decorate_plot(h,deco);
        
        subplot(1,2,2);
        % h = plot(Threshold_dB_MTC  ,Pcov_m_simul, 'r.',Threshold_dB_MTC  ,Pcov_m_analy, 'k-');
        hold on;
        for plt = 1:numel(LA_B)
            plot(params.Threshold.MTC_dB  ,Pcov_m_simul(plt,:), params.sim_style{plt},params.Threshold.MTC_dB  ,Pcov_m_analy(plt,:), params.ana_style{plt});
        end
        h = get(gca,'Children') ;
        
        deco.xlabel = 'Coverage threshold ($\tau$)  ';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = 'MTC';
        decorate_plot(h,deco);
        
        
        
        
        figure;plot(LA_B, IM,'r-',LA_B, IMH,'g--',LA_B,ISIC,'b:')
        figure;plot(LA_B, Bha,'r-',LA_B, Bhs,'ko');title('Bh');
        figure;plot(LA_B, Epha,'r-',LA_B, Ephs,'ko');title('Eph');
        figure;plot(LA_B, Epma,'r-',LA_B, Epms,'ko');title('Epm');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'Visualization';
        disp('Visualization');
        params.aggregation_mode = 'C2A' ;
        params.LA_B = 10000e-6;
        params.LA_H = 100e-6;
        
        
        visualize_the_network(params)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'Initial_Results';
        disp('Initial_Results');
        params.LA_M = 0.1:0.1:1;
        [R_oma , R_noma] = gen_initial_results(params);
        h = plot(params.LA_M *1e6 , R_oma * 1e-9 , 'k-' ,params.LA_M*1e6  ,R_noma* 1e-9  , 'r:');
        
        deco.xlabel = 'mMTC nodes density (nodes/km$^2$)  ';
        deco.ylabel = 'Network sum rate (Gbps/km$^2$) ';
        deco.title = '';
        deco.legend.items = [1 2] ;
        deco.legend.labels = {'OMA' , 'NOMA' };
        deco.legend.location = 'northwest';
        
        decorate_plot(h,deco);
end
end
