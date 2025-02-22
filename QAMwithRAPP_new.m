%% Simulation parameters 

EbN0dB = 0:1:20;          % Bit SNR range in dB
M16 = 16;                   % Example: 16-QAM
k16 = log2(M16);              % Bits per symbol
M64 = 64;
k64 = log2(M64);

PkNum=100; % initialise the number of packets to be transmitted


% Define RAPP PA model parameters
A0 = 1;     % Limiting output amplitude
v = 1;      % Small signal gain
p = 3;      % Smoothness parameter


%% Systems without AWGN 