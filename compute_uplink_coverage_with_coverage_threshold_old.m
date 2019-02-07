%function [Pcov_h, Pcov_m,params,Bha,Ep_h,Ep_m] = compute_uplink_coverage_with_coverage_threshold(params)
function [Pcov_h, Pcov_m,params] = compute_uplink_coverage_with_coverage_threshold_old(params)

th = params.Threshold.HTC;
tm = params.Threshold.MTC;
Pmo = params.MTC.Pmin;
Pmu = params.MTC.Pmax;
Pho = params.HTC.Pmin;
Phu = params.HTC.Pmax;
No =  params.No;
%theta_U = params.NOMA.theta_U;

 delta_h = th ./ (1 - th .* (Pmo / Pho)) ;
% delta_m = tm ./ (1 - tm .* (Pho / Pmo) .* theta_U) ;
 delta_m = tm ;


delta_h = delta_h(th < (Pho/Pmo));   % limit to the applicable range
%delta_m  = delta_m (tm < Pmo/(Pho * theta_U)  ); % limit to the applicable range
th = th(th < (Pho/Pmo));
%tm  = tm (tm < Pmo/(Pho * theta_U)  );
params.Threshold.HTC = th;
params.Threshold.HTC_dB = 10 * log10(th);
params.Threshold.MTC = tm;
params.Threshold.MTC_dB = 10 * log10(tm);


%  plot(th,delta_h)
%  figure;
%  plot(tm,delta_m)

a = params.SEPL.alpha;
b = params.SEPL.beta;
n = (2/b) - 1;

k = params.LA_B / params.LA_H;
po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
ph = 1 - po;
       
switch (params.aggregation_mode)
    case 'C2A'
        ps = ph;
    case 'C2C' 
        ps = 1;
end

theta_m = log((Pmu/Pmo).^(1/a)).^(1/b);
%ps = 1/2 * ps
Opower_m  = exp(- pi * ps .* params.LA_B .* theta_m^2);
theta_h = log((Phu/Pho).^(1/a)).^(1/b);
Opower_h  = exp(- pi .* params.LA_B .* theta_h^2);
%plot(params.LA_B * 1e6,Opower_m,params.LA_B * 1e6,Opower_h)

F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);

for k = 0:n
    NK = nchoosek(n,k);
    Eph(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta_h^2 , k/(n+1) + 1 , 'lower') ./ ((pi * params.LA_B) .^(k/(n+1)) .* (1 - Opower_h));
    Epm(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi * ps * params.LA_B .* theta_m^2 , k/(n+1) + 1 , 'lower') ./ ((pi * ps  .* params.LA_B) .^(k/(n+1)) .* (1 - Opower_m));
    for p = 1:numel(th)
        J_vh_h(k+1,p) = integral(@(y)F(y,n - k,delta_h(p)),0,delta_h(p));
        J_vh_m(k+1,p) = integral(@(y)F(y,n - k,delta_h(p)*(Pmo/Pho)),0,delta_h(p)*(Pmo/Pho));
    end
    for p = 1:numel(tm)
        J_vm_m(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)),0,delta_m(p));
        J_vm_h(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)*(Pho/Pmo)),0,delta_m(p)*(Pho/Pmo)); 
    end
    T_h_h(k+1,:) = NK * J_vh_h(k+1,:) * Eph(k+1) ;
    T_h_m(k+1,:) = NK * J_vh_m(k+1,:) * Epm(k+1) ;
    T_m_m(k+1,:) = NK * J_vm_m(k+1,:) * Epm(k+1) ;
    T_m_h(k+1,:) = NK * J_vm_h(k+1,:) * Eph(k+1) ;
end


E_h_h = sum(T_h_h,1);
E_h_m = sum(T_h_m,1);
E_m_m = sum(T_m_m,1);
E_m_h = sum(T_m_h,1);

A = (n+1) * pi / a^(n+1);
Bm =  params.rho_m * params.LA_M / params.N_RB  ;
Bh = ph*params.LA_B;

SNR_m = Pmo/No;
SNR_h = Pho/No;

noise_term_m = exp(- delta_m / SNR_m);
noise_term_h = exp(- delta_h / SNR_h);

interference_term_h_h = exp(- A * Bh * E_h_h);
interference_term_h_m = exp(- A * Bm * E_h_m);

interference_term_m_h = exp(- A * Bh * E_m_h);
interference_term_m_m = exp(- A * Bm * E_m_m);

Pcov_h =  noise_term_h .* interference_term_h_h .* interference_term_h_m;
Pcov_m =  noise_term_m .* interference_term_m_m .* interference_term_m_h; 

end