%% Scatter plot

constellation = '4PAM';
Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = pi/16;
frequency = 3*10^6;
hch = 1;

epsilon = 0.1;
hch_hat = epsilon*hch + hch;


starting_bits = round(rand(1,4000));
a = PHY.mapper(starting_bits, constellation);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);

demod = PHY.demodulate(mod, T, Ns, frequency, alpha, Lf, theta);

rdown = PHY.downsample(demod, Ns, Lf);

u = PHY.make_decision_variable(rdown,hch_hat,theta);

scatterplot(u);
title('');
hold on;

plot([0 0], [-1000 1000], 'Color', 'green')

hold off;

print('images/scatter_4PAM_aplitude','-depsc');
