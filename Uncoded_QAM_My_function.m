clc
clear
tic

%% initial parameters
M = 1000;          %number of symbols
mod_size = 16;     % modulation size, number of bits is log2(mod_size)
p = log2(mod_size);
info_bit = M*p;

SNR_min = 0;
SNR_max = 18;  % in dB
SNR_step = 1; % in dB
Run = 100;    % number of iteration in Monte-Carlo simulation

% Define RAPP PA model parameters
A0 = 1;     % Limiting output amplitude
v = 1;      % Small signal gain
prapp = 3;      % Smoothness parameter

%% Constellation generation
if mod_size == 4
    constellation = (1/sqrt(2)).*[-1-1j , -1+1j , 1-1j , 1+1j];  

    elseif mod_size == 16
    constellation = (1/sqrt(10)).*[-3 - 3j , -3 - 1j , -3 + 3j , -3 + 1j ,...
                     -1 - 3j , -1 - 1j , -1 + 3j , -1 + 1j ,...
                     3 - 3j , 3 - 1j , 3 + 3j , 3 + 1j ,...
                     1 - 3j , 1 - 1j , 1 + 3j , 1 + 1j];
else
    % This is for 64QAM but can be used to higher modulations like 256QAM
    % by changing the axis value
    axis = (1/sqrt(42)).*[-7 , -5 , -1 , -3 , 7 , 5 , 1 , 3];
    aaa = repmat(axis.' , [1 , sqrt(mod_size)]);
    bbb = repmat(axis , [sqrt(mod_size) , 1]);
    constellation = aaa + 1j*bbb;
    constellation = constellation.';
    constellation = constellation(:).';
end

%% constellation graph 
figure;
plot(real(constellation), imag(constellation), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('16-QAM Constellation');
%axis equal;


%% error vectors initialization
Error2 = zeros(((SNR_max - SNR_min)/SNR_step) + 1 , Run);

for run = 1 : Run
    run
%% Information bit generation
tx_info_bits = randi([0 1] , info_bit , 1); % generate random information bits

%% bit to symbol mapping
d1 = reshape(tx_info_bits , [p , length(tx_info_bits)/p]);         % produce matrix with p row at each column
decimal= bin2dec(num2str(d1.'));       % decimal value of bits
tx_symbols = constellation(decimal + 1).';        % symbols

%RAPP model addition 
modulated_signal_amp = abs(tx_symbols);
modulated_signal_phs = angle(tx_symbols);
modulated_signal_PA = RAPP_PA(modulated_signal_amp, A0, v, prapp);
modulated_signal_RAPP = modulated_signal_PA .* exp(1j * modulated_signal_phs); % Recombine

% tx_symbols = qammod(tx_info_bits , mod_size, "gray" , "InputType" , "bit" , "UnitAveragePower", true);

%% Adding AWGN noise 
Es = 1/p;    % average QAM symbol energy to bit energy

% Error vector initiations
Error1 = zeros(1 , ((SNR_max - SNR_min)/SNR_step) + 1);

cn = 1;   % counter for computing SNR
for SNR_dB = SNR_min : SNR_step : SNR_max
N0 = Es/(10^((SNR_dB)/10));
noise = sqrt(N0/2)*(randn(M , 1) + 1i*randn(M , 1));

%% Received signal
r = modulated_signal_RAPP + noise;

%% Symbol-to-bit mapping
   D = zeros(M , mod_size);
   % finding distance from all constellation points
   for ii = 1 : M
       for k = 1 : mod_size
           %D(ii , k) = distance([real(constellation(k)) , imag(constellation(k))] , [real(r(ii)) , imag(r(ii))]);
           D(ii, k) = abs(constellation(k) - r(ii));
       end
   end
   % Finding minimum index of distance e.g, 1 : p
   [~ , ind] = min(D ,[] , 2);

   bits_coded1 = (de2bi(ind - 1 , p,  'left-msb'))';    % changing numbers into binary. e.g. 3 = 11
   bits_hat = bits_coded1(:);

 % bits_hat = qamdemod(r , mod_size, "gray", "OutputType","bit", "UnitAveragePower", true);

%% Error rate calculation
Error1(cn) = (sum(tx_info_bits ~= bits_hat))/info_bit;

cn = cn + 1;
 end    % end for the SNR

 Error2(: , run) = Error1;

 end  % end for the run

%% Modulated signal 
figure;
plot(real(tx_symbols), imag(tx_symbols), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('16-QAM Signal');
SNR_plot = SNR_min : SNR_step : SNR_max;
figure;
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
% Add legend for clarity
legend('Theoretical BER','Simulated BER','Location','southwest')
toc
