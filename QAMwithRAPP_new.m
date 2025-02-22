%% Simulation parameters 

clc % clears the screen
clear all % clears all variables

EbN0dB = 0:1:20;          % Bit SNR range in dB
M16 = 16;                   % Example: 16-QAM
k16 = log2(M16);              % Bits per symbol
M64 = 64;
k64 = log2(M64);

numSymbols = 600; %number of symbols to be modulated 
PkNum=100; % initialise the number of packets to be transmitted


% Define RAPP PA model parameters
A0 = 1;     % Limiting output amplitude
v = 1;      % Small signal gain
p = 3;      % Smoothness parameter


%% Systems without AWGN 
% Data generetion 16 qam 
% Generate random binary data
dataBits16 = randi([0 1], numSymbols * k16, 1);
% Reshape the binary vector into a matrix where each row represents one symbol's bits
bit_matrix = reshape(dataBits16, k16, numSymbols).';
c = zeros(M16, 1);  % Pre-allocate the complex constellation array

% Convert each row (symbol) to a decimal index.
% For example, [0 0 1 1] becomes 3 (since '0011' in binary equals 3).
indices16 = bit_matrix * (2.^(k16-1:-1:0))';


%% 16 QAM model

for i = 0:M16-1

    % Convert i to a binary string of fixed length.
    % dec2bin returns a character array (e.g., '0011' for i=3 with num_bits_per_symbol=4)
    binStr16 = dec2bin(i, k16);
    % Convert the binary string to a numeric vector (0s and 1s)
    bits16 = binStr16 - '0';
        
    % For QAM, split the bits into two groups:
    %   bits(1:2:end) for the in-phase component
    %   bits(2:2:end) for the quadrature component
    pam_real = PAM_GRAY(bits16(1:2:end));
    pam_imag = PAM_GRAY(bits16(2:2:end));
        
    % Combine the two PAM components to form the QAM symbol
    c(i+1) = pam_real + 1i * pam_imag;
end

%Normalization
n = k16 / 2;
% Generate odd numbers: For n bits, these are 1, 3, ..., 2^n - 1.
odd_numbers = 1:2:(2^n - 1);
% Compute the variance factor (closed-form normalization factor):
qam_var = 1 / (2^(n-2)) * sum(odd_numbers.^2);
c16 = c / sqrt(qam_var);

%% Symbol mapping qam 16 
% Map each symbol index to its corresponding QAM constellation point.
% Since MATLAB indexing starts at 1, add 1 to the decimal index.
modulated_signal = c16(indices16 + 1);


%% RAPP PA addition 
modulated_signal_amp = abs(modulated_signal);
modulated_signal_phs = angle(modulated_signal);
modulated_signal_PA = RAPP_PA(modulated_signal_amp, A0, v, p);
modulated_signal_RAPP = modulated_signal_PA .* exp(1j * modulated_signal_phs); % Recombine
% Avoid division by zero: for nonzero elements, recombine amplitude and phase
nonzero = modulated_signal_amp > 0;
modulatedSignal_RAPP_2 = modulated_signal_PA(nonzero) .* (modulated_signal(nonzero)./modulated_signal_amp(nonzero));

%% Constellation plot 
figure;
plot(real(c16), imag(c16), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('16-QAM Constellation');
axis equal;

%% modulated signal plot 
figure;
plot(real(modulated_signal), imag(modulated_signal), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8); hold on;
plot(real(modulated_signal_RAPP), imag(modulated_signal_RAPP), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);
plot(real(modulatedSignal_RAPP_2), imag(modulatedSignal_RAPP_2), 'go', 'MarkerFaceColor','g', 'MarkerSize',8);
% grid on;
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('16-QAM Signal');
legend("Input Signal","RAPP PA Output","RAPP 2");
axis equal;



%% Demodulation: Minimum Distance Decision
% Pre-allocate vector for detected symbol indices
detected_indices = zeros(numSymbols,1);

% Loop over each received symbol
for n = 1:numSymbols
    % Compute squared Euclidean distances to all constellation points:
    distances = abs(modulated_signal(n) - c16).^2;
    % Find the index (MATLAB indexing: 1...M) of the minimum distance:
    [~, minIdx] = min(distances);
    % To be consistent with your modulation where indices started at 0,
    % we subtract 1:
    detected_indices(n) = minIdx - 1;
end

%% Map Detected Symbol Indices Back to Bits
% Pre-allocate bit vector for the demodulated bits:
demod_bits = zeros(numSymbols*k16, 1);

for n = 1:numSymbols
    % Convert the detected index to a binary string of length k
    binStr = dec2bin(detected_indices(n), k16);
    % Convert the string to a numeric vector (0s and 1s)
    bits = binStr - '0';
    % Place these bits in the appropriate location of the output vector
    demod_bits((n-1)*k16 + 1 : n*k16) = bits;
end

%% BER 
numErrors = sum(demod_bits ~= dataBits16);
BER = numErrors / length(dataBits16);
% %% RAPP signal plot 
% %% modulated signal plot 
% figure;
% plot(real(modulated_signal_RAPP), imag(modulated_signal_RAPP), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);
% grid on;
% xlabel('In-phase');
% ylabel('Quadrature');
% title('16-QAM RAPP Signal');
% axis equal;