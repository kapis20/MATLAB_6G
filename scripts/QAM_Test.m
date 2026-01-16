clc % clears the screen
clear all % clears all variables

M = 64;
k = log2(M);
%Create a binary data sequence. When using binary inputs, the number of rows in the input must be an integer multiple of the number of bits per symbol.

x = randi([0 1],10*k,1);
%Modulate the signal using bit inputs, and set it to have unit average power.

y = qammod(x,M,'bin', ...
    InputType='bit')%, ...
    %UnitAveragePower=true);
%Pass the signal through a noisy channel.
z = qamdemod(y,M,'bin',OutputType='bit');
s = isequal(x,double(z))
% rxSig = awgn(txSig,25);
%Plot the constellation diagram.

cd = comm.ConstellationDiagram(ShowReferenceConstellation=false);
cd(z)