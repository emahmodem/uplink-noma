        ind = 0;
        LA_B = [1000e-6 10000e-6 50000e-6];
        for i = 1:2:2*numel(LA_B)
             ind = ind + 1; 
             A{i} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
             A{i+1} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        
        A