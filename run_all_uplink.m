clc;close all;

subject = 'Simulation Progress';

scenarios = {   
'UL_Coverage_CoverageThreshold_C2A', ...
'UL_Coverage_SmallCellsDensity_C2A', ...
'UL_Coverage_HTC_Density_C2A', ...
'UL_Coverage_MTC_Density_C2A', ...
'UL_Coverage_NRB_C2A', ...
'UL_Coverage_Pho_Pmo_ratio_C2A', ...
'UL_Coverage_Pho_C2A'
};
for i = 1:length(scenarios)
    scenario = scenarios{i};
    disp(char(datetime))
    try
        gen_coverage_results_mMTC_NOMA_UDN_Uplink(scenario);
        send_report_via_email(subject,strcat(scenario ,' has been finished at: ' , char(datetime)));
    catch err
        disp(err)
        send_report_via_email(subject,strcat(scenario ,' has been failed at: ' , char(datetime),':',err.message));
    end  
end
 