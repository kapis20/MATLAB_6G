clc
clear
tic

%% initial parameters
M = 100;          %number of symbols
mod_size = 64;     % modulation size, number of bits is log2(mod_size)
p = log2(mod_size);
info_bit = M*p;

SNR_min = 0;
SNR_max = 18;  % in dB
SNR_step = 1; % in dB
Run = 100;    % number of iteration in Monte-Carlo simulation

%% error vectors initialization
Error2 = zeros(((SNR_max - SNR_min)/SNR_step) + 1 , Run);

for run = 1 : Run
    run
%% Information bit generation
tx_info_bits = randi([0 1] , info_bit , 1); % generate random information bits

tx_symbols = qammod(tx_info_bits , mod_size, "gray" , "InputType" , "bit" , "UnitAveragePower", true);

%% Adding AWGN noise 
Es = (1/log2(mod_size));    % average QAM symbol energy

% Error vector initiations
Error1 = zeros(1 , ((SNR_max - SNR_min)/SNR_step) + 1);

cn = 1;   % counter for computing SNR
for SNR_dB = SNR_min : SNR_step : SNR_max
N0 = Es/(10^((SNR_dB)/10));
noise = sqrt(N0/2)*(randn(M , 1) + 1i*randn(M , 1));

%% Received signal
r = tx_symbols + noise;

%% Symbol-to-bit mapping
 bits_hat = qamdemod(r , mod_size, "gray", "OutputType","bit", "UnitAveragePower", true);

%% Error rate calculation
Error1(cn) = (sum(tx_info_bits ~= bits_hat))/info_bit;

cn = cn + 1;
 end    % end for the SNR

 Error2(: , run) = Error1;

 end  % end for the run

SNR_plot = SNR_min : SNR_step : SNR_max;

%SER theory
[BER_theory , SER_theory] = berawgn(SNR_plot , 'QAM' , mod_size);
semilogy(SNR_plot , BER_theory)
hold on

%% BER
Error = mean(Error2 , 2);
semilogy(SNR_plot , Error , '--r')
hold on

xlabel('Eb/N0 (dB)')
ylabel('BER')
grid on

toc
