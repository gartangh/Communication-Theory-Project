function [ ] = scatter_plot(snr, constellation, m)
%SCATTER_PLOT Summary of this function goes here
%   Detailed explanation goes here
    Lf = 5;
    T = 1*10^(-6);
    Ns = 8;
    alpha = 1;
    theta = pi/16;
    frequency = 3*10^6;
    hch = 1;
    n = 1;
    k = 1;
    Eb = n/(k*m);
    N0 = Eb/snr;
    sigma = N0*Ns/(2*T);
    
    a = PHY.mapper(round(rand(1,1000)), constellation);
    
    mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
    signal = PHY.channel(mod, hch, sigma);


    demod = PHY.demodulate(signal, T, Ns, frequency, alpha, Lf, theta);

    rdown = PHY.downsample(demod, Ns, Lf);

    u = PHY.make_decision_variable(rdown,hch,theta);

    scatterplot(u);
    hold on;
    title(['Scatterplot voor E_B/N_0=' num2str(snr) 'db. Constellatie = ' constellation]);
    plot([0 0], [-1000 1000], 'Color', 'green')
    hold off;
    print([ 'images/scatter_' constellation '_' num2str(snr)],'-depsc');
end

