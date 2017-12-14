Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;

t = -Lf*T:T/Ns:Lf*T;
p = PHY.pulse(t,T,alpha);
plot(fourier(p));