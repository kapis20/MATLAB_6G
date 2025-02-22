% Complex Baseband BPSK Model
% AWGN Channel
% T. O'Farrell
% 2021.03.05

% Aim: to model BPSK in awgn

%initialisation
clc % clears the screen
clear all % clears all variables
randn('seed',0); % sets a seed for randn generator

% Define RAPP PA model parameters
A0 = 1;     % Limiting output amplitude
v = 1;      % Small signal gain
p = 3;      % Smoothness parameter

% Variable parameters
EbN0SdB=[0,1,2,3,4,5,6,7,8,9,10]; % initialise the EbN0 loop for simulation
EbN0TdB=[0,1,2,3,4,5,6,7,8,9,10]; % initialise the EbN0 loop for theory
PkLenBytes=100; % initialise the packet length in bytes
PkNum=1000; % initialise the number of packets to be transmitted
S=1; % initialise the transmit signal power

% Derived Parameters
PkLenBits=8*PkLenBytes; % data packet length in bits (convert from bytes to bits)
TxSignalLen=PkLenBits; % number of modulation symbols per packet

% EbN0 loop
for EbN0SIndex=1:length(EbN0SdB)

    % Loop derived parameters
    EbN0S=10^(EbN0SdB(EbN0SIndex)/10); % set EbN0 value for simulation
    StDev=sqrt(S/EbN0S); % set the noise standard deviation for calibration
    
% Pk Loop
    for PkIndex=1:PkNum
    
        % Transmitter
        TxData=rand(PkLenBits,1)>0.5; % generate the binary data
        TxSymbol=2*TxData-1; % BPSK data modulation
        
        %%%%%%%%%%%%%%% RAPP comes here %%%%%%%%%%%%%%%
        TxSymbol_amp = abs(TxSymbol);   % Get amplitudes
        TxSymbol_phase = angle(TxSymbol); % Get phases
        TxSymbol_modified_amp = RAPP_PA(TxSymbol_amp, A0, v, p); % Apply RAPP
        %Signal coming out of RAPP - multiply by the exponential form
        TxSymbol_RAPP = TxSymbol_modified_amp .* exp(1j * TxSymbol_phase); % Recombine
        % Complex baseband noise vector
        noise=StDev*(randn(TxSignalLen,1)+1i*randn(TxSignalLen,1))/sqrt(2);
        
        % Received signal vector
        h=1;
        %         h=(randn+1i*randn);
        RxSymbol=h*TxSymbol_RAPP + noise; % add noise to transmit signal
        
        % Receiver
        RxData=real(RxSymbol/h)>0; % zero threshold detection

        %Ber measurement
        %XOR compares transmitted and received data and sum adds up all the
        %errors 
        BitErrors(PkIndex)=sum(xor(TxData,RxData));
        PkErrors(PkIndex)=BitErrors(PkIndex)>0; %BLER
        
    end
    
    % Calculate error rates
    % XOR 
    ber(EbN0SIndex)=sum(BitErrors)/(PkNum*PkLenBits);
    per(EbN0SIndex)=sum(PkErrors)/PkNum;
end

% Theoretical BER and PER
for EbN0TIndex=1:length(EbN0TdB)
    EbN0T=10^(EbN0TdB(EbN0TIndex)/10); % set EbN0 value for theory
    tber(EbN0TIndex)=erfc(sqrt(EbN0T))/2; % theoretical BPSK BER performance
    tper(EbN0TIndex)=1-(1-tber(EbN0TIndex))^PkLenBits; % theoretical BPSK PER performance
end
    
% % Plot results
% figure
% semilogy(EbN0SdB,ber,'bd');
% hold
% semilogy(EbN0TdB,tber,'r-');
% figure
% semilogy(EbN0SdB,per,'bd');
% hold
% semilogy(EbN0TdB,tper,'r-');

% Plot results for BER
figure;
semilogy(EbN0SdB, ber, 'bd-', 'MarkerSize',8, 'LineWidth',1.5); hold on;
semilogy(EbN0TdB, tber, 'r-', 'LineWidth',1.5);
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER');
legend('Simulation', 'Theory');
grid on;

% Plot results for BLER (PER)
figure;
semilogy(EbN0SdB, per, 'bd-', 'MarkerSize',8, 'LineWidth',1.5); hold on;
semilogy(EbN0TdB, tper, 'r-', 'LineWidth',1.5);
xlabel('Eb/N0 (dB)');
ylabel('Block Error Rate (PER)');
title('BLER');
legend('Simulation', 'Theory');
grid on;

figure;
    
% Plot the waveform going into the RAPP PA (input signal)
subplot(2,1,1)
plot(real(TxSymbol), 'b.-');
title('Input Waveform to RAPP PA');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
    
% Plot the waveform coming out of the RAPP PA (output signal)
subplot(2,1,2)
plot(real(TxSymbol_RAPP), 'r.-');
title('Output Waveform from RAPP PA');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
% %AM/ AM 
% figure 
% plot(TxSymbol_amp,TxSymbol_modified_amp)
% xlabel('Input Amplitude');
% ylabel('Output Amplitude');
% title('RAPP PA AM/AM Conversion Curve');
% 
% figure
% pwelch(TxSymbol,[], [], [], 1, 'centered');
% title('PSD of the Signal using pwelch');
% xlabel('Normalized Frequency');
% ylabel('Power Spectral Density (dB/Hz)');



