function testDelta()
th_dB = (-20:0.1:20) ; % dB
th = 10.^(th_dB/10);
tm_dB = (-20:0.1:0) ; % dB
tm = 10.^(tm_dB/10);
e = 0.5;
theta_d = 0.05;

delta_h = th ./ (e - theta_d*(1-e).*th);
delta_m = tm ./ (1 - e.*(1 + tm));

figure;
plot(th_dB,delta_h,'k-')
figure;
plot(tm_dB,delta_m,'r.')

end