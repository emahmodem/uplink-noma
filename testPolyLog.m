function testPolyLog()
t_dB = (-20:0.5:20)  % dB
t = 10.^(t_dB/10)
n = 3
f1 = -abs(polylog(n,-t))
f2 = polylog(n,-t)
figure;
plot(t,f1,'k-',t,f2,'r.')

end