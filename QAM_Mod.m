%% initialisation
clc % clears the screen
clear all % clears all variables
randn('seed',0); % sets a seed for randn generator

% Simulation parameters
ebno_dB = 6:1:18;         % SNR range (in dB)
M = 64;                   % 64-QAM
bitsPerSymbol = log2(M);  % Should be 6 for 64-QAM
% numCodewords = 1000;      % Number of codewords (or packets)
% nBits = 4096;             % Codeword length in bits (adjust as needed)
% nSymbols = nBits / bitsPerSymbol; % Number of modulation symbols per codeword

PkLenBytesQAM = 510;           % New packet length in bytes -
PkNum = 1000;                 % Number of packets to be transmitted
PkLenBitsQAM = 8 * PkLenBytesQAM;  % 4080 (can divide by 6) and multiply of 8 - important! 
numSymbols = PkLenBitsQAM / bitsPerSymbol;     % Number of 64-QAM symbols per packet (4080/6 = 680)

% Pulse shaping filter parameters
rolloff = 0.3;
spanSymbols = 32;         % Filter span in symbols
sps = 4;                  % Samples per symbol (oversampling factor)

% RAPP PA parameters
A0 = 1;   % Saturation amplitude (same as before)
v  = 1;   % Small signal gain
p  = 3;   % Smoothness parameter

S=1; % initialise the transmit signal power

%% system 
%Noise loop 

for EbN0SIndex=1:length(ebno_dB)

    % Loop derived parameters
    EbN0S=10^(ebno_dB(EbN0SIndex)/10); % set EbN0 value for simulation
    StDev=sqrt(S/EbN0S); % set the noise standard deviation for calibration

    % packets loop:
    for PkIndex=1:PkNum
   
       % % Generate random data for the packet:
       % txBits = randi([0,M-1], PkLenBitsQAM, 1);
       %generate tandom bits 
       txBits = randi([0,1],PkLenBitsQAM,1);
      
       %Data modulation 
       txSymbol = qammod(txBits, M,'bin',InputType = 'bit',...
           UnitAveragePower=true);

        %%%%%%%%%%%%%%% RAPP comes here %%%%%%%%%%%%%%%
       
       % AWGN channel 

       % demodulation:
       rxbits = qamdemod(txSymbol,M,'bin',OutputType='bit');
       % s(EbN0SIndex) = isequal(txBits,double(rxbits))


       %Ber measurement
       %XOR compares transmitted and received data and sum adds up all the
       %errors 
       BitErrors(PkIndex)=sum(xor(txBits,rxbits));
       PkErrors(PkIndex)=BitErrors(PkIndex)>0; %BLER

       % Calculate error rates
       % XOR 
       ber(EbN0SIndex)=sum(BitErrors)/(PkNum*PkLenBitsQAM);
       per(EbN0SIndex)=sum(PkErrors)/PkNum;
    end

end

scatterplot(txSymbol)

% Plot results for BER
figure;
semilogy(ebno_dB, ber, 'bd', 'MarkerSize',8, 'LineWidth',1.5); hold on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER');
legend('Simulation')
grid on;

% Plot results for BLER (PER)
figure;
semilogy(ebno_dB, per, 'bd', 'MarkerSize',8, 'LineWidth',1.5); hold on;
xlabel('Eb/N0 (dB)');
ylabel('Block Error Rate (PER)');
title('BLER');
legend('Simulation')
grid on;