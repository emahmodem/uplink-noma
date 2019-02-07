Pmo = 1e-11;
Pmu = 1e-3;
Pho = 1e-10;
Phu =10e-3;
No =  1e-13;
a = 0.3;
b = 2/3;
n = (2/b) - 1;
LA_B = [1e-6 10e-6 100e-6 1000e-6 10000e-6 50000e-6 100000e-6];
LA_H = 500e-6;
k = LA_B / LA_H;
po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
ph = 1 - po;
aggregation_mode =  'C2C';      
switch (aggregation_mode)
    case 'C2A'
        ps = ph;
    case 'C2C' 
        ps = 1;
end



%params.LA_B = [0.001e-6:0.001e-6:1e-5];
theta_m = log((Pmu/Pmo).^(1/a)).^(1/b);
%ps = 1/2 * ps
Opower_m  = exp(- pi * ps .* LA_B .* theta_m^2);
theta_h = log((Phu/Pho).^(1/a)).^(1/b);
Opower_h  = exp(- pi .* LA_B .* theta_h^2);

k = 1;

    Eph =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* LA_B .* theta_h^2 , k/(n+1) + 1 , 'lower') ./ ((pi * LA_B) .^(k/(n+1)) .* (1 - Opower_h));
    Epm =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi * ps .* LA_B .* theta_m^2 , k/(n+1) + 1 , 'lower') ./ ((pi * ps  .*LA_B) .^(k/(n+1)) .* (1 - Opower_m));
    
Bh = ph.*LA_B;

figure;plot(LA_B * 1e6,Opower_m,'-r',LA_B * 1e6,Opower_h,'--b');
figure;plot(LA_B * 1e6,Epm,'-r',LA_B * 1e6,Eph,'--b');
figure;plot(LA_B * 1e6,Bh,'--b');

% F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% 
% 
% 
% 
% 
% for k = 0:n
%     NK = nchoosek(n,k);
%     Eph(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta_h^2 , k/(n+1) + 1 , 'lower') ./ ((pi * params.LA_B) .^(k/(n+1)) .* (1 - Opower_h));
%     Epm(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi * ps * params.LA_B .* theta_m^2 , k/(n+1) + 1 , 'lower') ./ ((pi * ps  .* params.LA_B) .^(k/(n+1)) .* (1 - Opower_m));
%     for p = 1:numel(th)
%         J_vh_h(k+1,p) = integral(@(y)F(y,n - k,delta_h(p)),0,delta_h(p));
%         J_vh_m(k+1,p) = integral(@(y)F(y,n - k,delta_h(p)),0,delta_h(p)*(Pmo/Pho));
%     end
%     for p = 1:numel(tm)
%         J_vm_m(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)),0,delta_m(p));
%         J_vm_h(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)),0,delta_m(p)*(Pho/Pmo));
%     end
%     T_h_h(k+1,:) = NK * J_vh_h(k+1,:) * Eph(k+1) ;
%     T_h_m(k+1,:) = NK * J_vh_m(k+1,:) * Eph(k+1) ;
%     T_m_m(k+1,:) = NK * J_vm_m(k+1,:) * Epm(k+1) ;
%     T_m_h(k+1,:) = NK * J_vm_h(k+1,:) * Epm(k+1) ;
% end
% 
% E_h_h = sum(T_h_h,1);
% E_h_m = sum(T_h_m,1);
% E_m_m = sum(T_m_m,1);
% E_m_h = sum(T_m_h,1);
% 
% A = (n+1) * pi / a^(n+1);
% Bm =  params.rho_m * params.LA_M / params.N_RB  ;
% Bh = ph*params.LA_B;
% 
% SNR_m = Pmo/No;
% SNR_h = Pho/No;
% 
% noise_term_m = exp(- delta_m / SNR_m);
% noise_term_h = exp(- delta_h / SNR_h);
% 
% interference_term_h_h = exp(- A * Bh * E_h_h);
% interference_term_h_m = exp(- A * Bm * E_h_m);
% interference_term_m_h = exp(- A * Bh * E_m_h);
% interference_term_m_m = exp(- A * Bm * E_m_m);
% 
% Pcov_h = noise_term_h .* interference_term_h_h .* interference_term_h_m;
% Pcov_m = noise_term_m .* interference_term_m_m .* interference_term_m_h; 
% 
