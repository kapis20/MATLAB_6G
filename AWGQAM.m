%% Simulation parameters 

clc % clears the screen
clear all % clears all variables

EbN0SdB = 0:1:20;          % Bit SNR range in dB
M16 = 16;                   % Example: 16-QAM
k16 = log2(M16);              % Bits per symbol
M64 = 64;
k64 = log2(M64);

numSymbols = 600; %number of symbols to be modulated 
PkNum=1000; % initialise the number of packets to be transmitted

S=1; % initialise the transmit signal power

% Define RAPP PA model parameters
A0 = 1;     % Limiting output amplitude
v = 1;      % Small signal gain
p = 3;      % Smoothness parameter


%% AWGN simualtion loop 
Es = 1/k16;
% EbN0 loop
for EbN0SIndex=1:length(EbN0SdB)

    % Loop derived parameters
    EbN0S=10^(EbN0SdB(EbN0SIndex)/10); % set EbN0 value for simulation
    StDev=sqrt(S/EbN0S); % set the noise standard deviation for calibration
    N0 = Es/(10^((EbN0SIndex)/10));
    % Calculate noise variance per dimension for a complex baseband signal
    % Eb/N0 (linear) is EbN0S, k is bits per symbol, and Es=1 is assumed.
    sigma = sqrt(1/(2*k16*EbN0S));  % standard deviation per real dimension

    for PkIndex=1:PkNum
        % Generate random binary data
        dataBits = randi([0 1], numSymbols * k16, 1);
        % Reshape the binary vector into a matrix where each row represents one symbol's bits
        bit_matrix = reshape(dataBits, k16, numSymbols).';
        c = zeros(M16, 1);  % Pre-allocate the complex constellation array
        
        % Convert each row (symbol) to a decimal index.
        % For example, [0 0 1 1] becomes 3 (since '0011' in binary equals 3).
        indices = bit_matrix * (2.^(k16-1:-1:0))';
        %%% QAM Model %%% 
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

        % constellation normalization 
        c16 = c/sqrt(10);

        %% Symbol mapping qam 16 
        % Map each symbol index to its corresponding QAM constellation point.
        % Since MATLAB indexing starts at 1, add 1 to the decimal index.
        modulated_signal = c16(indices + 1);

        %% RAPP PA addition 
        modulated_signal_amp = abs(modulated_signal);
        modulated_signal_phs = angle(modulated_signal);
        modulated_signal_PA = RAPP_PA(modulated_signal_amp, A0, v, p);
        modulated_signal_RAPP = modulated_signal_PA .* exp(1j * modulated_signal_phs); % Recombine

        %% AWGN channel 
        % Complex baseband noise vector
        noise1=StDev*(randn(numSymbols,1)+1i*randn(numSymbols,1))/sqrt(2);
        
        % Generate complex AWGN noise vector
        noise = sigma * (randn(numSymbols,1) + 1i*randn(numSymbols,1));
        %noise = sqrt(N0/2)*(randn(numSymbols , 1) + 1i*randn(numSymbols , 1));
        % Received signal vector
        h=1;
        %         h=(randn+1i*randn);
        RxSymbols=h*modulated_signal+ noise; % add noise to transmit signal
        
        RxSymbolsRAPP = h*modulated_signal_RAPP+ noise;

        %% Demodulation: Minimum Distance Decision
        % Pre-allocate vector for detected symbol indices
        detected_indices = zeros(numSymbols,1);

        % Loop over each received symbol
        for n = 1:numSymbols
            % Compute squared Euclidean distances to all constellation points:
            distances = abs(RxSymbols(n) - c16).^2;
            % Find the index (MATLAB indexing: 1...M) of the minimum distance:
            [~, minIdx] = min(distances);
            % To be consistent with your modulation where indices started at 0,
            % we subtract 1:
            detected_indices(n) = minIdx - 1;
        end

        detected_indicesRapp=zeros(numSymbols,1);
        for n = 1:numSymbols
            % Compute squared Euclidean distances to all constellation points:
            distancesRAPP = abs(RxSymbolsRAPP(n) - c16).^2;
            % Find the index (MATLAB indexing: 1...M) of the minimum distance:
            [~, minIdx] = min(distancesRAPP);
            % To be consistent with your modulation where indices started at 0,
            % we subtract 1:
            detected_indicesRapp(n) = minIdx - 1;
        end

        %% Map Detected Symbol Indices Back to Bits
        % Pre-allocate bit vector for the demodulated bits:
        demod_bits = zeros(numSymbols*k16, 1);
        demod_bitsRapp = zeros(numSymbols*k16, 1);
        for n = 1:numSymbols
            % Convert the detected index to a binary string of length k
            binStr = dec2bin(detected_indices(n), k16);
            % Convert the string to a numeric vector (0s and 1s)
            bits = binStr - '0';
            % Place these bits in the appropriate location of the output vector
            demod_bits((n-1)*k16 + 1 : n*k16) = bits;
        end

        detected_indicesRapp=zeros(numSymbols,1);
        for n = 1:numSymbols
            % Convert the detected index to a binary string of length k
            binStr = dec2bin(detected_indicesRapp(n), k16);
            % Convert the string to a numeric vector (0s and 1s)
            bits = binStr - '0';
            % Place these bits in the appropriate location of the output vector
            demod_bitsRapp((n-1)*k16 + 1 : n*k16) = bits;
        end


   
     
        BERerrors(PkIndex) = sum(demod_bits ~= dataBits);
        %BERerrors(PkIndex) = numErrors / length(dataBits);
        BLERErrors(PkIndex)=BERerrors(PkIndex)>0; %BLER
        
        BERerrorsRapp(PkIndex) = sum(demod_bitsRapp ~= dataBits);
        %BERerrors(PkIndex) = numErrors / length(dataBits);
        BLERErrorsRapp(PkIndex)=BERerrorsRapp(PkIndex)>0; %BLER
    end

    BER(EbN0SIndex)=sum(BERerrors)/ (PkNum*length(dataBits));
    BLER(EbN0SIndex)=sum(BLERErrors)/PkNum;

    BERRapp(EbN0SIndex)=sum(BERerrorsRapp)/ (PkNum*length(dataBits));
    BLERRapp(EbN0SIndex)=sum(BLERErrorsRapp)/PkNum;
end


%% theorhethical AWGN QAM 
% Convert Eb/N0 from dB to linear scale
EbN0lin = 10.^(EbN0SdB/10);
%k bit average energy 
EsN0lin16 = k16 .* EbN0lin;  % gamma_s = k * gamma_b
Ps16 = zeros(size(EbN0SdB));

% Compute the symbol error probability using the erfc version of the formula
for i = 1:length(EbN0SdB)
    %gamma_s256 = EsN0lin256(i);
    gamma_s16 = EsN0lin16(i);
    Ps16(i) = 2*erfc(sqrt(3*gamma_s16/(2*(M16-1))))
    %Ps64(i) = 2*erfc(sqrt(3*gamma_s64/(2*(M64-1))))
end

% Approximate BER from SER using Gray coding approximation
BER16 = Ps16 / k16;

% Compute Block Error Rate (BLER)
BLER16 = 1 - (1 - BER16).^PkNum;
%BLER64 = 1 - (1 - BER64).^block_size;
%BER64 = Ps64 / k64;

%% BER, BLER plots 
figure;
semilogy(EbN0SdB, BER, 'b-o', 'LineWidth', 2); hold on;
semilogy(EbN0SdB, BERRapp, 'g-^', 'LineWidth', 2); 
semilogy(EbN0SdB, BER16, 'r-s', 'LineWidth', 2); 
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
ylim([1e-5 1e-1]);  % Adjust the y-axis limits as needed
legend('BER', 'BER RAPP', 'Theoretical BER');
title('16-QAM BER Performance');

figure;
semilogy(EbN0SdB, BLER, 'b-s', 'LineWidth', 2); hold on;
semilogy(EbN0SdB, BLERRapp, 'm-d', 'LineWidth', 2);
semilogy(EbN0SdB, BLER16, 'r-s', 'LineWidth', 2); 
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BLER');
ylim([1e-5 1e-1]);  % Adjust the y-axis limits as needed
legend('BLER', 'BLER RAPP','Theoretical BLER');
title('16-QAM BLER Performance');

