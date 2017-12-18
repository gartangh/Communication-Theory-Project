
BER = [];
correct_BER = [];
starting_bits = round(rand(1,300000));
a = PHY.mapper(starting_bits, constellation);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
for snr = 1:0.1:8
    constellation = '4QAM';
    Lf = 5;
    T = 1*10^(-6);
    Ns = 8;
    alpha = 1;
    theta = pi/16;
    frequency = 3*10^6;
    hch = 1;
    n = 1;
    k = 1;
    m = 2;
    Eb = n/(k*m);
    N0 = Eb/snr;
    sigma = N0*Ns/(2*T);

    signal = PHY.channel(mod, hch, sigma);


    demod = PHY.demodulate(signal, T, Ns, frequency, alpha, Lf, theta);

    rdown = PHY.downsample(demod, Ns, Lf);

    u = PHY.make_decision_variable(rdown,hch,theta);

    estim = PHY.hard_decisions(u, constellation);

    result_bits = PHY.demapper(estim, constellation);

    BER = [BER (sum(result_bits~=starting_bits)/length(result_bits))];
    correct_BER = [correct_BER qfunc(sqrt(2*snr))];
    disp(snr);
end

semilogy(mag2db(1:0.1:8)/2, BER);
hold on;
semilogy(mag2db(1:0.1:8)/2, correct_BER);

xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('Simulatie', 'Analytisch');
hold off;
print('images/BER_4QAM','-depsc');
