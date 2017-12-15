Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = 0;
frequency = 3*10^(-6);
hch = 1;
constellation = 'BPSK';
a = PHY.mapper(round(rand(1,100)), constellation);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
signal = PHY.channel(mod, hch, theta);
demod = PHY.demodulate(mod, T, Ns, frequency, alpha, Lf, theta);
result = demod(2*Lf*Ns+1:length(demod) - 2*Lf*Ns);
plot(result);
rdown = PHY.downsample(demod, Ns, Lf);

u = PHY.make_decision_variable(rdown,hch,theta);

estim = PHY.hard_decisions(u, constellation);
disp(abs(estim-a)<1e-6);