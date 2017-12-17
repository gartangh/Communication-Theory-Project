%% Constants
Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = pi/16;
frequency = 3*10^6;
hch = 1;
 
% Constellation constants
constellation = 'BPSK';
m = 2;
n = 1;
k = 1;
snr = 5.4;

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


%disp(q);

%% Scenario 1
%samples_bits2 = Channel_Coding.Encode_outer(samples_bits);
%samples_bits2 = samples_bits;
%% Scenario 2

a = PHY.mapper(samples_bits, constellation);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
signal = PHY.channel(mod, hch, sigma);


demod = PHY.demodulate(signal, T, Ns, frequency, alpha, Lf, theta);

rdown = PHY.downsample(demod, Ns, Lf);

u = PHY.make_decision_variable(rdown,hch,0);

estim = PHY.hard_decisions(u, constellation);
%disp(length(a)-sum(a==estim));


result_bits = PHY.demapper(estim, constellation);


%result_bits2 = Channel_Coding.Decode_outer(result_bits);

result = zeros(1, length(result_bits)/NBITS);
for i = 0:length(result)-1
    result(i+1) = bi2de(result_bits(i*NBITS+1:(i+1)*NBITS))+1;
end
result = arrayfun(@(i) q(i), result);
Quantization.show_figures(result);

% Aantal verschillende bits
disp(sum(result~=samples_quantized));