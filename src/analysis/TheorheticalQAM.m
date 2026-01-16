% Define simulation parameters
EbN0dB = 0:1:20;          % Bit SNR range in dB
M256 = 256;
k256 = log2(M256)
M16 = 16;                   % Example: 16-QAM
k16 = log2(M16);              % Bits per symbol
M64 = 64;
k64 = log2(M64);


% Convert Eb/N0 from dB to linear scale
EbN0lin = 10.^(EbN0dB/10);
% 
% Compute symbol SNR: Es/N0 = k * Eb/N0
%EsN0lin256 = k256 .* EbN0lin;  % gamma_s = k * gamma_b
%k bit average energy 
EsN0lin16 = k16 .* EbN0lin;  % gamma_s = k * gamma_b
EsN0lin64 = k64 .* EbN0lin;  % gamma_s = k * gamma_b
% Preallocate array for symbol error probability
Ps256 = zeros(size(EbN0dB));
Ps16 = zeros(size(EbN0dB));
Ps64 = zeros(size(EbN0dB));
% Compute the symbol error probability using the erfc version of the formula
for i = 1:length(EbN0dB)
    %gamma_s256 = EsN0lin256(i);
    gamma_s16 = EsN0lin16(i);
    gamma_s64 = EsN0lin64(i);
    % 
    % term1= 2*(1-1/sqrt(M));
    % term2 = term1*erfc(sqrt(3* gamma_s/(2*(M-1))));
    % term3 = term2* (1-1/2*(1-1/sqrt(M)));
    % 
    % Ps(i) = term3*erfc(sqrt(3*gamma_s/(2*(M-1))));

    %Ps4(i) = 2*erfc(sqrt(3*gamma_s256/(2*(M-1))))
    Ps16(i) = 2*erfc(sqrt(3*gamma_s16/(2*(M16-1))))
    Ps64(i) = 2*erfc(sqrt(3*gamma_s64/(2*(M64-1))))
end

% Plot the symbol error probability versus Eb/N0
figure;
%semilogy(EbN0dB, Ps256,  'r-o', 'LineWidth', 2); hold on;
semilogy(EbN0dB, Ps16, 'b-s', 'LineWidth', 2); hold on;
semilogy(EbN0dB, Ps64, 'g-d', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Symbol error');
% Set the y-axis limit to only show values up to 10^-5
ylim([1e-5 1e-1]);  % Adjust the lower limit as needed
title('Theoretical SER for 16-, and 64-QAM using erfc');
legend('16-QAM', '64-QAM', 'Location', 'southwest');


SNR_min = 1;
SNR_step = 1;
SNR_max = 20;


SNR_plot = SNR_min : SNR_step : SNR_max;

[BER_theory16 , SER_theory] = berawgn(SNR_plot , 'QAM' , M16);
[BER_theory64 , SER_theory] = berawgn(SNR_plot , 'QAM' , M64);

% Approximate BER from SER using Gray coding approximation
BER16 = Ps16 / k16;
BER64 = Ps64 / k64;
figure;
semilogy(EbN0dB, BER16, 'b-s', 'LineWidth', 2); hold on;
semilogy(EbN0dB, BER64, 'g-d', 'LineWidth', 2);
semilogy(SNR_plot, BER_theory16, 'r-s', 'LineWidth', 2);
semilogy(SNR_plot, BER_theory64, 'y-d', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BLER)');
ylim([1e-5 1e-1]);  % Adjust the y-axis limits as needed
title('Theoretical Bit Error Rate for 16-QAM and 64-QAM');
legend('16-QAM', '64-QAM','16- QAM MAT','64 MAT', 'Location', 'southwest');



% Define block size (in bits)
block_size = 600;

% Compute Block Error Rate (BLER)
BLER16 = 1 - (1 - BER16).^block_size;
BLER64 = 1 - (1 - BER64).^block_size;
% Plot the block error rates versus Eb/N0
figure;
semilogy(EbN0dB, BLER16, 'b-s', 'LineWidth', 2); hold on;
semilogy(EbN0dB, BLER64, 'g-d', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');
ylim([1e-5 1e-1]);  % Adjust the y-axis limits as needed
title('Theoretical Block Error Rate for 16-QAM and 64-QAM');
legend('16-QAM', '64-QAM', 'Location', 'southwest');