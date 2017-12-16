Ns = 8;
T = 10^-6;
frequency = 3*10^6;
alpha = 1;
Lf = 5;
theta = pi/16;
h_ch = 1;
s2n = 5.4;
h_ch_hat = 1;
theta_hat = pi/16;

t_afgeknot = -Lf*T:T/Ns:Lf*T;
s_afgeknot = PHY.pulse(t_afgeknot, T, alpha);

t_normaal = -10*Lf*T:T/Ns:10*Lf*T;
s_normaal = PHY.pulse(t_normaal, T, alpha);


fourier_afgeknot = fft(s_afgeknot);
frequentie_afgeknot = (0:length(fourier_afgeknot)-1)/T/length(fourier_afgeknot);

fourier_normaal = fft(s_normaal);
frequentie_normaal = (0:length(fourier_normaal)-1)/T/length(fourier_normaal);

plot(frequentie_afgeknot,abs(fourier_afgeknot));
hold on;
plot(frequentie_normaal,abs(fourier_normaal));
hold off;
xlabel('Frequentie (Hz)');
ylabel('DFT (genormaliseerd)');

%% Nu nog een focking uitleg :(