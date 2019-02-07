th= 0.0001:0.0001:9;
Pmo = 1;
Pho = 10;

delta_h = th ./ (1 - th .* (Pmo / Pho)) ;


plot(th,delta_h)