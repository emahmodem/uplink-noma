%function [Pcov_h, Pcov_m,params,Bha,Ep_h,Ep_m] = compute_uplink_coverage_with_coverage_threshold(params)
function [Pcov_h, Pcov_m] = compute_uplink_coverage_with_coverage_threshold(params)

th = params.Threshold.HTC;
tm = params.Threshold.MTC;
th_qos = params.Threshold.HTC_QOS;
Pmo = params.MTC.Pmin;
Pmu = params.MTC.Pmax;
Pho = params.HTC.Pmin;
Phu = params.HTC.Pmax;
No =  params.No;
%theta_U = params.NOMA.theta_U;

 delta_h = th ./ (1 - th .* (Pmo / Pho)) ;
% delta_m = tm ./ (1 - tm .* (Pho / Pmo) .* theta_U) ;
 delta_m = tm ;

delta_h_qos = th_qos ./ (1 - th_qos .* (Pmo / Pho)) ;
epsi = th_qos + Pho/Pmo * tm + th_qos * tm;
%delta_h = delta_h(th < (Pho/Pmo));   % limit to the applicable range
%delta_m  = delta_m (tm < Pmo/(Pho * theta_U)  ); % limit to the applicable range
%th = th(th < (Pho/Pmo));
%tm  = tm (tm < Pmo/(Pho * theta_U)  );
% params.Threshold.HTC = th;
% params.Threshold.HTC_dB = 10 * log10(th);
% params.Threshold.MTC = tm;
% params.Threshold.MTC_dB = 10 * log10(tm);


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
    
    
    [~, po] = find(th < Pho/Pmo,1,'last') ;
    
    for p = 1:po
        J_vh_h_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h(p)),0,2*delta_h(p));
        J_vh_h_r1_t2(k+1,p) = integral(@(y)F(y,n - k,th(p)),0,th(p));
        J_vh_m_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h(p)*(Pmo/Pho)),0,2*delta_h(p)*(Pmo/Pho));
        J_vh_m_r1_t2(k+1,p) = integral(@(y)F(y,n - k,th(p)*(Pmo/Pho)),0,th(p)*(Pmo/Pho));
    end
    
    for p = (po+1) :numel(th)
        J_vh_h_r2(k+1,p-po) = integral(@(y)F(y,n - k,th(p)),0,th(p));
        J_vh_m_r2(k+1,p-po) = integral(@(y)F(y,n - k,th(p)*(Pmo/Pho)),0,th(p)*(Pmo/Pho));
    end
        

    if(th_qos * Pmo/Pho - 1.0 >=0 )
        for p = 1:numel(tm)
            J_vm_m_T1(k+1,p) = integral(@(y)F(y,n - k,epsi(p)*(Pmo/Pho)),0,epsi(p)*(Pmo/Pho));
            J_vm_h_T1(k+1,p) = integral(@(y)F(y,n - k,epsi(p)),0,epsi(p));
        end
    else
        J_vm_m_T1(k+1,p) = 0;
        J_vm_h_T1(k+1,p) = 0;
        [~, pom] = find(tm < Pmo/Pho * delta_h_qos,1,'last') ;
        for p = 1:pom
            J_vm_m_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h_qos*(Pmo/Pho)),0,2*delta_h_qos*(Pmo/Pho));
            J_vm_m_r1_t2(k+1,p) = integral(@(y)F(y,n - k,epsi(p)*(Pmo/Pho)),0,epsi(p)*(Pmo/Pho));
            J_vm_h_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h_qos),0,2*delta_h_qos);
            J_vm_h_r1_t2(k+1,p) = integral(@(y)F(y,n - k,epsi(p)),0,epsi(p));
        end
        
        for p = (pom+1) :numel(tm)
            J_vm_m_r2(k+1,p-pom) = integral(@(y)F(y,n - k,2*tm(p)),0,2*tm(p));
            J_vm_h_r2(k+1,p-pom) = integral(@(y)F(y,n - k,2*tm(p)*(Pho/Pmo)),0,2*tm(p)*(Pho/Pmo));
        end
    end
%     for p = 1:numel(tm)
%         J_vm_m(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)),0,delta_m(p));
%         J_vm_h(k+1,p) = integral(@(y)F(y,n - k,delta_m(p)*(Pho/Pmo)),0,delta_m(p)*(Pho/Pmo)); 
%     end
    
 

    
    T_h_h_r1_t1(k+1,:) = NK * J_vh_h_r1_t1(k+1,:) * Eph(k+1) ;
    T_h_h_r1_t2(k+1,:) = NK * J_vh_h_r1_t2(k+1,:) * Eph(k+1) ;
    T_h_h_r2(k+1,:) = NK * J_vh_h_r2(k+1,:) * Eph(k+1) ;
    T_h_m_r1_t1(k+1,:) = NK * J_vh_m_r1_t1(k+1,:) * Epm(k+1) ;
    T_h_m_r1_t2(k+1,:) = NK * J_vh_m_r1_t2(k+1,:) * Epm(k+1) ;
    T_h_m_r2(k+1,:) = NK * J_vh_m_r2(k+1,:) * Epm(k+1) ;
    
    T_m_m_T1(k+1,:) = NK * J_vm_m_T1(k+1,:) * Epm(k+1) ;
    T_m_h_T1(k+1,:) = NK * J_vm_h_T1(k+1,:) * Eph(k+1) ;
    
    T_m_m_r1_t1(k+1,:) = NK * J_vm_m_r1_t1(k+1,:) * Epm(k+1) ;
    T_m_m_r1_t2(k+1,:) = NK * J_vm_m_r1_t2(k+1,:) * Epm(k+1) ;
    T_m_m_r2(k+1,:) = NK * J_vm_m_r2(k+1,:) * Epm(k+1) ;
    T_m_h_r1_t1(k+1,:) = NK * J_vm_h_r1_t1(k+1,:) * Eph(k+1) ;
    T_m_h_r1_t2(k+1,:) = NK * J_vm_h_r1_t2(k+1,:) * Eph(k+1) ;
    T_m_h_r2(k+1,:) = NK * J_vm_h_r2(k+1,:) * Eph(k+1) ;
end


E_h_h_r1_t1 = sum(T_h_h_r1_t1,1);
E_h_h_r1_t2 = sum(T_h_h_r1_t2,1);
E_h_h_r2 = sum(T_h_h_r2,1);
E_h_m_r1_t1 = sum(T_h_m_r1_t1,1);
E_h_m_r1_t2 = sum(T_h_m_r1_t2,1);
E_h_m_r2 = sum(T_h_m_r2,1);

E_m_m_T1 = sum(T_m_m_T1,1);
E_m_h_T1 = sum(T_m_h_T1,1);

E_m_m_r1_t1 = sum(T_m_m_r1_t1,1);
E_m_m_r1_t2 = sum(T_m_m_r1_t2,1);
E_m_m_r2 = sum(T_m_m_r2,1);
E_m_h_r1_t1 = sum(T_m_h_r1_t1,1);
E_m_h_r1_t2 = sum(T_m_h_r1_t2,1);
E_m_h_r2 = sum(T_m_h_r2,1);

A = (n+1) * pi / a^(n+1);
Bm =  params.rho_m * params.LA_M / params.N_RB  ;
Bh = ph*params.LA_B;

SNR_m = Pmo/No;
SNR_h = Pho/No;



noise_term_h_r1_t1 = exp(-2* delta_h(1:po) / SNR_h);
noise_term_h_r1_t2 = exp(-th(1:po) / SNR_h);
noise_term_h_r2 = exp(-th(po+1:end) / SNR_h);


interference_term_h_h_r1_t1 = exp(- A * Bh * E_h_h_r1_t1);
interference_term_h_h_r1_t2 = exp(- A * Bh * E_h_h_r1_t2);
interference_term_h_h_r2 = exp(- A * Bh * E_h_h_r2);
interference_term_h_m_r1_t1 = exp(- A * Bm * E_h_m_r1_t1);
interference_term_h_m_r1_t2 = exp(- A * Bm * E_h_m_r1_t2);
interference_term_h_m_r2 = exp(- A * Bm * E_h_m_r2);

noise_term_m_T1 = exp(- epsi / SNR_h);
noise_term_m_r1_t1 = exp(-2* delta_h(1:pom) / SNR_h);
noise_term_m_r1_t2 = exp(-epsi(1:pom) / SNR_h);
noise_term_m_r2 = exp(-2*tm(pom+1:end) / SNR_m);

interference_term_m_h_T1 = exp(- A * Bh * E_m_h_T1);
interference_term_m_m_T1 = exp(- A * Bm * E_m_m_T1);

interference_term_m_m_r1_t1 = exp(- A * Bm * E_m_m_r1_t1);
interference_term_m_m_r1_t2 = exp(- A * Bm * E_m_m_r1_t2);
interference_term_m_m_r2 = exp(- A * Bm * E_m_m_r2);
interference_term_m_h_r1_t1 = exp(- A * Bh * E_m_h_r1_t1);
interference_term_m_h_r1_t2 = exp(- A * Bh * E_m_h_r1_t2);
interference_term_m_h_r2 = exp(- A * Bh * E_m_h_r2);

r = th*Pmo/Pho;
T1 = (r-1)./(r+1);
T2 = 2./(1+r);

r_qos = th_qos*Pmo/Pho;
T1m = (r_qos-1)./(r_qos+1);
T2m = 2./(1+r_qos);

Pcov_h_r1_t1 =  T1(1:po) .* noise_term_h_r1_t1 .* interference_term_h_h_r1_t1 .* interference_term_h_m_r1_t1;
Pcov_h_r1_t2 =  T2(1:po).* noise_term_h_r1_t2 .* interference_term_h_h_r1_t2 .* interference_term_h_m_r1_t2;
Pcov_h_r2 =  T2(po+1:end) .* noise_term_h_r2 .* interference_term_h_h_r2 .* interference_term_h_m_r2;

Pcov_m_r1_t1 =  T1m * noise_term_m_r1_t1 .* interference_term_m_m_r1_t1 .* interference_term_m_h_r1_t1;
Pcov_m_r1_t2 =  T2m * noise_term_m_r1_t2 .* interference_term_m_m_r1_t2 .* interference_term_m_h_r1_t2;
Pcov_m_r2 =   noise_term_m_r2 .* interference_term_m_m_r2 .* interference_term_m_h_r2;


Pcov_h = [Pcov_h_r1_t1 + Pcov_h_r1_t2  Pcov_h_r2];

if(th_qos * Pmo/Pho - 1.0 >=0)
    Pcov_m =  T2m*noise_term_m_T1 .* interference_term_m_m_T1 .* interference_term_m_h_T1; 
else
    Pcov_m = [Pcov_m_r1_t1 + Pcov_m_r1_t2  Pcov_m_r2];
end
       
end