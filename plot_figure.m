close all;
%load DL_Coverage_CoverageThreshold_C2A ; fig = 'DL_Coverage_CoverageThreshold_C2A';
%load DL_Coverage_SmallCellsDensity_C2A; fig = 'DL_Coverage_SmallCellsDensity_C2A';
load DL_Coverage_NOMA_PC_C2A; fig = 'DL_Coverage_NOMA_PC_C2A';
switch(fig)
    case 'DL_Coverage_NOMA_PC_C2A'
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

        savefig('Downlink_Initial_Results/dl_coverage_epsilon_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_epsilon_c2a','-depsc','-r0');
        
    case 'DL_Coverage_SmallCellsDensity_C2A'
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
        
        savefig('Downlink_Initial_Results/dl_coverage_las_c2c');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_las_c2c','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_las_c2c','-depsc','-r0');
        
    case 'DL_Coverage_CoverageThreshold_C2A'
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
        savefig('Downlink_Initial_Results/dl_coverage_tau_c2a');
        set(gcf,'PaperPositionMode','auto')
        print('Downlink_Initial_Results/dl_coverage_tau_c2a','-dpng','-r0');
        print('Downlink_Initial_Results/dl_coverage_tau_c2a','-depsc','-r0');
end