function [Pcov_h, Pcov_m] = compute_uplink_coverage_with_Pho_Pmo_ratio(params)
th = params.Threshold.HTC;
tm = params.Threshold.MTC;
th_qos = params.Threshold.HTC_QOS;
Pmo = params.MTC.Pmin;
Pmu = params.MTC.Pmax;
Pho = params.HTC.Pmin;
Phu = params.HTC.Pmax;
No =  params.No;

delta_h = th ./ (1 - th .* (Pmo ./ Pho)) ;
delta_h_qos = th_qos ./ (1 - th_qos .* (Pmo ./ Pho)) ;
epsi = th_qos + Pho./Pmo * tm + th_qos * tm;


a = params.SEPL.alpha;
b = params.SEPL.beta;
n = (2/b) - 1;

k = params.LA_B ./ params.LA_H;
po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
ph = 1 - po;

switch (params.aggregation_mode)
    case 'C2A'
        ps = ph;
    case 'C2C'
        ps = 1;
end

theta_m = log((Pmu./Pmo).^(1/a)).^(1/b);
%ps = 1/2 * ps
Opower_m  = exp(- pi * ps .* params.LA_B .* theta_m.^2);
theta_h = log((Phu./Pho).^(1/a)).^(1/b);
Opower_h  = exp(- pi .* params.LA_B .* theta_h.^2);
%plot(params.LA_B * 1e6,Opower_m,params.LA_B * 1e6,Opower_h)

F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);

for k = 0:n
    NK = nchoosek(n,k);
    Eph(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta_h^2 , k/(n+1) + 1 , 'lower') ./ ((pi * params.LA_B) .^(k/(n+1)) .* (1 - Opower_h));
    Epm(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi * ps .* params.LA_B .* theta_m.^2 , k/(n+1) + 1 , 'lower') ./ ((pi * ps  .* params.LA_B) .^(k/(n+1)) .* (1 - Opower_m));
    for p = 1:numel(Pmo)
        if(th < Pho/Pmo(p))
            J_vh_h_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h(p)),0,2*delta_h(p));
            J_vh_h_r1_t2(k+1,p) = integral(@(y)F(y,n - k,th),0,th);
            J_vh_m_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h(p).*(Pmo(p)/Pho)),0,2*delta_h(p)*(Pmo(p)/Pho));
            J_vh_m_r1_t2(k+1,p) = integral(@(y)F(y,n - k,th*(Pmo(p)/Pho)),0,th*(Pmo(p)/Pho));
            
            T_h_h_r1_t1(k+1,p) = NK * J_vh_h_r1_t1(k+1,p) * Eph(k+1) ;
            T_h_h_r1_t2(k+1,p) = NK * J_vh_h_r1_t2(k+1,p) * Eph(k+1) ;
            T_h_m_r1_t1(k+1,p) = NK * J_vh_m_r1_t1(k+1,p) * Epm(k+1,p) ;
            T_h_m_r1_t2(k+1,p) = NK * J_vh_m_r1_t2(k+1,p) * Epm(k+1,p) ;
        else
            J_vh_h_r2(k+1,p) = integral(@(y)F(y,n - k,th),0,th);
            J_vh_m_r2(k+1,p) = integral(@(y)F(y,n - k,th*(Pmo(p)/Pho)),0,th*(Pmo(p)/Pho));
            
            T_h_h_r2(k+1,p) = NK * J_vh_h_r2(k+1,p) * Eph(k+1) ;
            T_h_m_r2(k+1,p) = NK * J_vh_m_r2(k+1,p) * Epm(k+1,p) ;
        end
        
        if(th_qos * Pmo(p)/Pho - 1.0 >=0 )
            J_vm_m_T1(k+1,p) = integral(@(y)F(y,n - k,epsi(p)*(Pmo(p)/Pho)),0,epsi(p)*(Pmo(p)/Pho));
            J_vm_h_T1(k+1,p) = integral(@(y)F(y,n - k,epsi(p)),0,epsi(p));
            
            T_m_m_T1(k+1,p) = NK * J_vm_m_T1(k+1,p) * Epm(k+1,p) ;
            T_m_h_T1(k+1,p) = NK * J_vm_h_T1(k+1,p) * Eph(k+1) ;
        else
            if(tm < Pmo(p)/Pho * delta_h_qos(p))
                J_vm_m_T1(k+1,p) = 0;
                J_vm_h_T1(k+1,p) = 0;
                J_vm_m_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h_qos(p)*(Pmo(p)/Pho)),0,2*delta_h_qos(p)*(Pmo(p)/Pho));
                J_vm_m_r1_t2(k+1,p) = integral(@(y)F(y,n - k,epsi(p)*(Pmo(p)/Pho)),0,epsi(p)*(Pmo(p)/Pho));
                J_vm_h_r1_t1(k+1,p) = integral(@(y)F(y,n - k,2*delta_h_qos(p)),0,2*delta_h_qos(p));
                J_vm_h_r1_t2(k+1,p) = integral(@(y)F(y,n - k,epsi(p)),0,epsi(p));
                
                T_m_m_r1_t1(k+1,p) = NK * J_vm_m_r1_t1(k+1,p) * Epm(k+1,p) ;
                T_m_m_r1_t2(k+1,p) = NK * J_vm_m_r1_t2(k+1,p) * Epm(k+1,p) ;
                T_m_h_r1_t1(k+1,p) = NK * J_vm_h_r1_t1(k+1,p) * Eph(k+1) ;
                T_m_h_r1_t2(k+1,p) = NK * J_vm_h_r1_t2(k+1,p) * Eph(k+1) ;
            else
                J_vm_m_r2(k+1,p) = integral(@(y)F(y,n - k,2*tm),0,2*tm);
                J_vm_h_r2(k+1,p) = integral(@(y)F(y,n - k,2*tm*(Pho/Pmo(p))),0,2*tm*(Pho/Pmo(p)));
                
                T_m_m_r2(k+1,p) = NK * J_vm_m_r2(k+1,p) * Epm(k+1,p) ;
                T_m_h_r2(k+1,p) = NK * J_vm_h_r2(k+1,p) * Eph(k+1) ;
            end
        end
    end
    %     T_h_h_r1_t1(k+1,:) = NK * J_vh_h_r1_t1(k+1,:) * Eph(k+1) ;
    %     T_h_h_r1_t2(k+1,:) = NK * J_vh_h_r1_t2(k+1,:) * Eph(k+1) ;
    %     T_h_h_r2(k+1,:) = NK * J_vh_h_r2(k+1,:) * Eph(k+1) ;
    %     T_h_m_r1_t1(k+1,:) = NK * J_vh_m_r1_t1(k+1,:) * Epm(k+1) ;
    %     T_h_m_r1_t2(k+1,:) = NK * J_vh_m_r1_t2(k+1,:) * Epm(k+1) ;
    %     T_h_m_r2(k+1,:) = NK * J_vh_m_r2(k+1,:) * Epm(k+1) ;
    
    %     T_m_m_T1(k+1,:) = NK * J_vm_m_T1(k+1,:) * Epm(k+1) ;
    %     T_m_h_T1(k+1,:) = NK * J_vm_h_T1(k+1,:) * Eph(k+1) ;
    %
    %     T_m_m_r1_t1(k+1,:) = NK * J_vm_m_r1_t1(k+1,:) * Epm(k+1) ;
    %     T_m_m_r1_t2(k+1,:) = NK * J_vm_m_r1_t2(k+1,:) * Epm(k+1) ;
    %     T_m_m_r2(k+1,:) = NK * J_vm_m_r2(k+1,:) * Epm(k+1) ;
    %     T_m_h_r1_t1(k+1,:) = NK * J_vm_h_r1_t1(k+1,:) * Eph(k+1) ;
    %     T_m_h_r1_t2(k+1,:) = NK * J_vm_h_r1_t2(k+1,:) * Eph(k+1) ;
    %     T_m_h_r2(k+1,:) = NK * J_vm_h_r2(k+1,:) * Eph(k+1) ;
end

A = (n+1) * pi / a^(n+1);
Bm =  params.rho_m .* params.LA_M ./ params.N_RB  ;
Bh = ph.*params.LA_B;

SNR_m = Pmo./No;
SNR_h = Pho./No;

noise_term_h_r1_t1 = exp(-2* delta_h ./ SNR_h);
noise_term_h_r1_t2 = exp(-th ./ SNR_h);
noise_term_h_r2 = exp(-th ./ SNR_h);

noise_term_m_T1 = exp(- epsi ./ SNR_h);
noise_term_m_r1_t1 = exp(-2* delta_h ./ SNR_h);
noise_term_m_r1_t2 = exp(-epsi ./ SNR_h);
noise_term_m_r2 = exp(-2*tm ./ SNR_m);

r = th.*Pmo./Pho;
T1 = (r-1)./(r+1);
T2 = 2./(1+r);

r_qos = th_qos.*Pmo./Pho;
T1m = (r_qos-1)./(r_qos+1);
T2m = 2./(1+r_qos);

for p = 1:numel(Pmo)
    if(th < Pho/Pmo(p))
        E_h_h_r1_t1(p) = sum(T_h_h_r1_t1(:,p));
        E_h_h_r1_t2(p) = sum(T_h_h_r1_t2(:,p));
        E_h_m_r1_t1(p) = sum(T_h_m_r1_t1(:,p));
        E_h_m_r1_t2(p) = sum(T_h_m_r1_t2(:,p));
        interference_term_h_h_r1_t1(p) = exp(- A * Bh .* E_h_h_r1_t1(p));
        interference_term_h_h_r1_t2(p) = exp(- A * Bh .* E_h_h_r1_t2(p));
        interference_term_h_m_r1_t1(p) = exp(- A * Bm .* E_h_m_r1_t1(p));
        interference_term_h_m_r1_t2(p) = exp(- A * Bm .* E_h_m_r1_t2(p));
        Pcov_h_r1_t1(p) =  T1(p) .* noise_term_h_r1_t1(p) .* interference_term_h_h_r1_t1(p) .* interference_term_h_m_r1_t1(p);
        Pcov_h_r1_t2(p) =  T2(p).* noise_term_h_r1_t2 .* interference_term_h_h_r1_t2(p) .* interference_term_h_m_r1_t2(p);
        Pcov_h = Pcov_h_r1_t1 + Pcov_h_r1_t2  ;
    else
        E_h_h_r2(p) = sum(T_h_h_r2(:,p));
        E_h_m_r2(p) = sum(T_h_m_r2(:,p));
        interference_term_h_h_r2(p) = exp(- A * Bh .* E_h_h_r2(p));
        interference_term_h_m_r2(p) = exp(- A * Bm .* E_h_m_r2(p));
        Pcov_h_r2(p) =  T2(p) .* noise_term_h_r2(p) .* interference_term_h_h_r2(p) .* interference_term_h_m_r2(p);
        Pcov_h(p) = Pcov_h_r2(p);
    end
    
    if(th_qos * Pmo(p)/Pho - 1.0 >=0 )
        E_m_m_T1(p) = sum(T_m_m_T1,1);
        E_m_h_T1(p) = sum(T_m_h_T1,1);
        interference_term_m_h_T1(p) = exp(- A * Bh .* E_m_h_T1(p));
        interference_term_m_m_T1(p) = exp(- A * Bm .* E_m_m_T1(p));
        Pcov_m(p) =  T2m(p)*noise_term_m_T1(p) .* interference_term_m_m_T1(p) .* interference_term_m_h_T1(p);
    else
        if(tm < Pmo(p)/Pho * delta_h_qos(p))
            E_m_m_r1_t1(p) = sum(T_m_m_r1_t1(:,p));
            E_m_m_r1_t2(p) = sum(T_m_m_r1_t2(:,p));
            E_m_h_r1_t1(p) = sum(T_m_h_r1_t1(:,p));
            E_m_h_r1_t2(p) = sum(T_m_h_r1_t2(:,p));
            interference_term_m_m_r1_t1(p) = exp(- A * Bm .* E_m_m_r1_t1(p));
            interference_term_m_m_r1_t2(p) = exp(- A * Bm .* E_m_m_r1_t2(p));
            interference_term_m_h_r1_t1(p) = exp(- A * Bh .* E_m_h_r1_t1(p));
            interference_term_m_h_r1_t2(p) = exp(- A * Bh .* E_m_h_r1_t2(p));
            Pcov_m_r1_t1(p) =  T1m(p) * noise_term_m_r1_t1(p) .* interference_term_m_m_r1_t1(p) .* interference_term_m_h_r1_t1(p);
            Pcov_m_r1_t2(p)=  T2m(p) * noise_term_m_r1_t2(p) .* interference_term_m_m_r1_t2(p) .* interference_term_m_h_r1_t2(p);
            Pcov_m(p) = Pcov_m_r1_t1(p) + Pcov_m_r1_t2(p) ;
        else
            E_m_m_r2(p) = sum(T_m_m_r2(:,p));
            E_m_h_r2(p) = sum(T_m_h_r2(:,p));
            interference_term_m_m_r2(p) = exp(- A * Bm .* E_m_m_r2(p));
            interference_term_m_h_r2(p) = exp(- A * Bh .* E_m_h_r2(p));
            Pcov_m_r2(p) =   noise_term_m_r2(p) .* interference_term_m_m_r2(p) .* interference_term_m_h_r2(p);
            Pcov_m(p) =  Pcov_m_r2(p);
        end
    end
end

end