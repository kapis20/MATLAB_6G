% Define simulation parameters
EbN0dB = 0:1:20;                % Bit SNR range in dB
M_values = [4, 16, 64];         % Modulation orders for QAM
colors = ['r', 'b', 'g'];       % Colors for each curve
markers = ['o', 's', 'd'];      % Markers for each curve

figure; hold on;

for idx = 1:length(M_values)
    M = M_values(idx);
    k = log2(M);               % Bits per symbol

    % Convert Eb/N0 from dB to linear scale
    EbN0lin = 10.^(EbN0dB/10);

    % Compute symbol SNR: Es/N0 = k * Eb/N0
    EsN0lin = k .* EbN0lin;     % gamma_s = k * gamma_b

    % Preallocate array for symbol error probability
    Ps = zeros(size(EbN0dB));

    % Compute the symbol error probability using the erfc version of the formula
    for i = 1:length(EbN0dB)
        gamma_s = EsN0lin(i);
        Ps(i) = 2 * erfc( sqrt(3 * gamma_s / (2*(M-1))) );
    end

    % Plot the SER curve for this modulation order
    semilogy(EbN0dB, Ps, [colors(idx) '-' markers(idx)], 'LineWidth', 2);
end

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Symbol Error Probability (SER)');
ylim([1e-5 1e-1]);  % Set the y-axis limits
title('Theoretical SER for 4-, 16-, and 64-QAM using erfc');
legend('4-QAM', '16-QAM', '64-QAM','Location','southwest');
