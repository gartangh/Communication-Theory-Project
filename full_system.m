%% Constants
Lf = 5;
T = 1*10^(-6);
Ns = 8;
alpha = 1;
theta = pi/16;
frequency = 3*10^(-6);
hch = 1;
constellation = 'BPSK';

%% Kwantisatie

% Get Lloyd Max data
%[GKD,SQR,entropie,r,q,p] = Quantization.Lloyd_max_quantizer;
% Get linear quantizer data
[Delta_opt,GKD,SQR,entropie,r,q,p] = Quantization.optimal_linear_quantizer;


% Quantize the figure
[samples_quantized]= Quantization.quantize(r,q);
% Show side by side comparison
%Quantization.show_figures(samples_quantized);


%{

MAX van de samples = 188 met 1 na de komma => 1880 => 11 bits per 

%}


NBITS = 11;
samples_bits = zeros(1,NBITS*length(samples_quantized));
for i=0:length(samples_quantized)-1
    samples_bits(i*NBITS+1:(i+1)*NBITS) = de2bi(uint16(samples_quantized(i+1)*10), 11);
end


%% Scenario 1
samples_bits = Channel_Coding.Encode_outer(samples_bits);

%% Scenario 2
disp('a');

a = PHY.mapper(samples_bits, constellation);

mod = PHY.modulate(a,T,Ns,frequency,alpha,Lf);
signal = PHY.channel(mod, hch, theta);
demod = PHY.demodulate(mod, T, Ns, frequency, alpha, Lf, theta);
result = demod(2*Lf*Ns+1:length(demod) - 2*Lf*Ns);
%plot(result);
rdown = PHY.downsample(demod, Ns, Lf);

u = PHY.make_decision_variable(rdown,hch,theta);

estim = PHY.hard_decisions(u, constellation);
disp('b');
%disp(length(a)-sum(a==estim));


result_bits = PHY.demapper(estim, constellation);
result = zeros(1, length(result_bits)/NBITS);


for i = 0:length(result)-1
    result(i+1) = bi2de(result_bits(i*NBITS+1:(i+1)*NBITS))/10.0;
end

result = Channel_Coding.Decode_outer(result);

Quantization.show_figures(result);