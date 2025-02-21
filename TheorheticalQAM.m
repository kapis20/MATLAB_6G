% Define simulation parameters
EbN0dB = 0:1:20;          % Bit SNR range in dB
M = 16;                   % Example: 16-QAM
k = log2(M);              % Bits per symbol


% Compute symbol SNR: Es/N0 = k * Eb/N0
%gamma = EbN0dB .* k;

%convert to linear 
%gamma = 10.^(gamma/10);
% Convert Eb/N0 from dB to linear scale
EbN0lin = 10.^(EbN0dB/10);
% 
% Compute symbol SNR: Es/N0 = k * Eb/N0
EsN0lin = k .* EbN0lin;  % gamma_s = k * gamma_b

% Preallocate array for symbol error probability
Ps = zeros(size(EbN0dB));

% Compute the symbol error probability using the erfc version of the formula
for i = 1:length(EbN0dB)
    gamma_s = EsN0lin(i);
    % 
    % term1= 2*(1-1/sqrt(M));
    % term2 = term1*erfc(sqrt(3* gamma_s/(2*(M-1))));
    % term3 = term2* (1-1/2*(1-1/sqrt(M)));
    % 
    % Ps(i) = term3*erfc(sqrt(3*gamma_s/(2*(M-1))));

    Ps(i) = 2*erfc(sqrt(3*gamma_s/(2*(M-1))))
end

% Plot the symbol error probability versus Eb/N0
figure;
semilogy(EbN0dB, Ps, 'r-o', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Symbol Error Probability (SER)');
% Set the y-axis limit to only show values up to 10^-5
ylim([1e-5 1e-1]);  % Adjust the lower limit as needed
title(sprintf('Theoretical SER for %d-QAM using erfc', M));
