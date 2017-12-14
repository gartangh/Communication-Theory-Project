Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = pi/16;
frequency = 3*10^(-6);

a = PHY.mapper(round(rand(1,1000)), 'BPSK');
disp(a);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
demod = PHY.demodulate(mod, T, Ns, frequency, alpha, Lf, theta);
result = demod(2*Lf*Ns+1:length(demod) - 2*Lf*Ns);

eyediagram(real(result), 2*Ns, 2*Ns)