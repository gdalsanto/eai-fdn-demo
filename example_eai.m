close all; clear all; clc

addpath(genpath('utility'))
addpath(genpath('fdnToolbox'))

fs = 441000;    % sampling frequency
N = 8;          % matrix size
irLen = fs*2;   % length of impulse response 
t = (0:1:irLen-1)/fs;       % time samples

A = randomOrthogonal(N);    % orthogonal matrix 
B = ones(N,1);              % input gains
C = ones(1,N);              % output gains
D = zeros(1,1);             % direct gain
delays = [809, 877, 937, 1049, 1151, 1249, 1373, 1499]; % delays

% RT = 0.5;                   % target reverberation time
% g = 10^(-3/fs/RT);          % attenuation per sample
g = 1;
Gamma = diag(g.^delays);    
A = A*Gamma;                % apply homogeneous decay

%% EAI

% use EAI method to compute poles 
[residues, poles, direct, isConjugatePolePair, metaData] = ...
    dss2pr(delays,A, B, C, D);  

% generate impulse response
ir = dss2impz(irLen, delays, A, B, C, D);

%% plot
% impulse response
figure('Name','impulse response');
plot(t, ir);
xlabel('time [s]');
ylabel('amplitude');
grid on 
title('FDN impulse response - lossless')

% poles
figure('Name','poles');
plot(real(poles),imag(poles),'x');
xlabel('real');
ylabel('imag');
grid on 
title('FDN system poles - lossless')

%% apply absorption filter
% absorption filters
centerFrequencies = [ 63, 125, 250, 500, 1000, 2000, 4000, 8000]; % Hz
T60frequency = [1, centerFrequencies fs];
targetT60 = [4; 4; 4.2; 4.3; 5; 6; 6; 6; 8; 10];  % seconds
% absorption filter 
zAbsorption = zSOS(absorptionGEQ(targetT60, delays, fs),'isDiagonal',true);

% power correction filter
targetPower = [5; 5; 5; 3; 2; 1; -1; -3; -5; -5];  % dB
powerCorrectionSOS = designGEQ(targetPower);
C = zSOS(permute(powerCorrectionSOS,[3 4 1 2]) .* C);

% use EAI method to compute poles 
[residues, poles, direct, isConjugatePolePair, metaData] = ...
    dss2pr(delays,A, B, C, D,'absorptionFilters', zAbsorption);  

% generate impulse response
ir = dss2impz(irLen, delays, A, B, C, D,'absorptionFilters', zAbsorption);

rt_poles = -60/fs/20./log10(abs(poles));

%% plot
% impulse response
figure('Name','impulse response with absorption');
plot(t, ir);
xlabel('time [s]');
ylabel('amplitude');
grid on 
title('FDN impulse response - with absorption')

% poles
norm_rt_poles = rt_poles/max(rt_poles);
figure('Name','poles with absorption');
plot(norm_rt_poles.*real(poles./abs(poles)),norm_rt_poles.*imag(poles./abs(poles)),'x');
hold on
r = linspace(0, pi, 100);
plot(cos(r), sin(r),'.')
xlabel('real');
ylabel('imag');
grid on 
title('FDN system poles (normalized RT)')