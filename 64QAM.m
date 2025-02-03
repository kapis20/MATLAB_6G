%% initialisation
clc % clears the screen
clear all % clears all variables
randn('seed',0); % sets a seed for randn generator

% Simulation parameters
ebno_dB = 6:1:18;         % SNR range (in dB)
% M = 64;                   % 64-QAM
% bitsPerSymbol = log2(M);  % Should be 6 for 64-QAM
% numCodewords = 1000;      % Number of codewords (or packets)
% nBits = 4096;             % Codeword length in bits (adjust as needed)
% nSymbols = nBits / bitsPerSymbol; % Number of modulation symbols per codeword

PkLenBytesQAM = 96;           % New packet length in bytes -
PkNum = 1000;                 % Number of packets to be transmitted
PkLenBitsQAM = 8 * PkLenBytes;  % 768 bits, and 768/6 = 128 symbols exactly
numSymbols = PkLenBits / 6;     % Number of 64-QAM symbols per packet (768/6 = 128)

% Pulse shaping filter parameters
rolloff = 0.3;
spanSymbols = 32;         % Filter span in symbols
sps = 4;                  % Samples per symbol (oversampling factor)

% RAPP PA parameters
A0 = 1;   % Saturation amplitude (same as before)
v  = 1;   % Small signal gain
p  = 3;   % Smoothness parameter

S=1; % initialise the transmit signal power

%% 