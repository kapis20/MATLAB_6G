%PSD Calculations 
fft_size = 2* PkLenBits;
% Compute the FFT and center the zero-frequency component
fft_TxData = fft(TxData, fft_size);
fft_TxData = fftshift(fft_TxData);


% Calculate the Power Spectral Density (PSD)
psd = (abs(fft_TxData).^2) / fft_size;  % Normalized by the signal length

% Create a normalized frequency axis (range -0.5 to 0.5)
f_axis = linspace(-0.5, 0.5, fft_size);

% Plot the PSD in dB scale
figure;
plot(f_axis, 10*log10(psd), 'LineWidth', 1.5);
xlabel('Normalized Frequency');
ylabel('Power (dB)');
title('Power Spectral Density of TxData (FFT Method)');
grid on;