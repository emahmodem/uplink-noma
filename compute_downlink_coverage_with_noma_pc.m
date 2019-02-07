function [Pcov_h, Pcov_m] = compute_downlink_coverage_with_noma_pc(params)

a = params.SEPL.alpha;
b = params.SEPL.beta;
th = params.Threshold.HTC;
tm = params.Threshold.MTC;
n = (2/b) - 1;
e = params.NOMA.epsilon;
theta_d = params.NOMA.theta_D;
delta_h = th ./ (e - theta_d.*(1-e).*th);
delta_m = tm ./ (1 - e.*(1 + tm));
NN =  params.No / params.PD ;

k = 0:n;
T = (pi * factorial(n+1)) ./ (factorial(k) .* a.^(n-k+1));
for k = 0:n
    apl_h(k+1,:) = T(k+1) .* (- abs(polylog(n-k+1,-delta_h)));
    apl_m(k+1,:) = T(k+1) .* (- abs(polylog(n-k+1,-delta_m)));
end

switch(b)
    case 1/2
        F = @(r,ps,pa,la_b,delta,NN,a,b,apl,ind)  r .* exp(- pi * ps .* la_b .* r.^2  ) .* exp (- delta .* NN .* exp(a .* r .^ b))  ...
            .* exp(pa * la_b .* apl(1,ind)) .* exp(pa * la_b .* apl(2,ind).* r.^(0.5)) ...
            .* exp(pa * la_b .* apl(3,ind).* r) .* exp(pa * la_b .* apl(4,ind).* r.^(1.5));
    case 2/3
        F = @(r,ps,pa,la_b,delta,NN,a,b,apl,ind)  r .* exp(- pi * ps .* la_b .* r.^2  ).* exp (- delta .* NN .* exp(a .* r .^ b))  ...
            .* exp(pa * la_b .* apl(1,ind)) .* exp(pa * la_b .* apl(2,ind).* r.^(2/3)) ...
            .* exp(pa * la_b .* apl(3,ind).* r.^(4/3));
    case 1
        F = @(r,ps,pa,la_b,delta,NN,a,b,apl,ind)  r .* exp(- pi * ps .* la_b .* r.^2  ).* exp (- delta .* NN .* exp(a .* r .^ b))  ...
            .* exp(pa * la_b .* apl(1,ind)) .* exp(pa * la_b .* apl(2,ind).* r) ;
        
end

switch (params.aggregation_mode)
    case 'C2A'
        k = params.LA_B ./ params.LA_H;
        po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
        pa = 1 - po;
        ps = pa;
    case 'C2C' 
        k= params.LA_B ./ (params.LA_H + params.rho_m * params.LA_M);
        po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
        pa = 1 - po;
        ps = 1;
end

if(b == 2)
    Pcov_h = exp(- pi * pa .* (params.LA_B ./ a) .* log(1 + delta_h));
    Pcov_m = exp(- pi * pa .* (params.LA_B ./ a) .* log(1 + delta_m));
    return;
end


for p = 1:numel(params.NOMA.epsilon)
    Pcov_h(p) =  2 * pi * params.LA_B * integral(@(r)F(r,1,pa,params.LA_B,delta_h(p),NN,a,b,apl_h,p),0,inf);
    Pcov_m(p) =  2 * pi *  ps *  params.LA_B * integral(@(r)F(r,ps,pa,params.LA_B,delta_m(p),NN,a,b,apl_m,p),0,inf);
end

end