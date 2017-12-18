%% Constants

Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = pi/16;
frequency = 3*10^6;
hch = 1;
constellation = 'BPSK';
m = 1;
n = 1;
k = 1;
snr = 2.5;

%% Calculate sigma

Eb = n/(k*m);
N0 = Eb/snr;
sigma = N0*Ns/(2*T);
disp(sigma);

%% Kwantisatie

% Get Lloyd Max data
[GKD,SQR,entropie,r,q,p] = Quantization.Lloyd_max_quantizer;
% Get linear quantizer data
%[Delta_opt,GKD,SQR,entropie,r,q,p] = Quantization.optimal_linear_quantizer;


% Quantize the figure
[samples_quantized]= Quantization.quantize(r,q);

samples_quantized_idx = arrayfun(@(x)find(q==x,1),samples_quantized);

NBITS = 3;
samples_bits = zeros(1,NBITS*length(samples_quantized_idx));
for i=0:length(samples_quantized_idx)-1
    samples_bits(i*NBITS+1:(i+1)*NBITS) = de2bi(uint16(samples_quantized_idx(i+1))-1, NBITS);
end

result_bits_decoded = zeros(1, length(samples_bits));
counts = zeros(1, length(samples_bits)/5);
for k = 0:length(samples_bits)/5-1
    bool_error = 1;
    while bool_error
        counts(k+1) = counts(k+1)+1;
        samples_bits_part = Channel_Coding.Encode_inner(samples_bits(k*5+1:(k+1)*5), [1 1 0 1 0 1]);

        a = PHY.mapper(samples_bits_part, constellation);

        mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
        signal = PHY.channel(mod, hch, sigma);


        demod = PHY.demodulate(signal, T, Ns, frequency, alpha, Lf, theta);

        rdown = PHY.downsample(demod, Ns, Lf);

        u = PHY.make_decision_variable(rdown,hch,theta);

        estim = PHY.hard_decisions(u, constellation);

        result_bits = PHY.demapper(estim, constellation);
        [result_bits_decoded_part, bool_error] = Channel_Coding.Decode_inner(result_bits, [1 1 0 1 0 1]);
        result_bits_decoded(k*5+1:(k+1)*5) = result_bits_decoded_part;
    end
end
result = zeros(1, length(result_bits_decoded)/NBITS);
for i = 0:length(result)-1
    result(i+1) = bi2de(result_bits_decoded(i*NBITS+1:(i+1)*NBITS))+1;
end
result = arrayfun(@(i) q(i), result);
Quantization.show_figures(result);
% Aantal verschillende bits
BER = sum(result~=samples_quantized)/length(result);

[original, ~] = make_probability_functions(Quantization.filename);
original = reshape(original, 1, length(result));
GKA = mean((result-original).^2);



disp(['BER = ' num2str(BER)]);
disp(['Gemiddeld aantal bits = ' num2str(sum(counts)*10/length(result))]);
disp(['GKA = ' num2str(GKA)]);

%{
BER = 0
Gemiddeld aantal bits = 7.1075
GKA = 35.5276
%}