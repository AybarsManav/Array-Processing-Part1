function [X,A,S] = gendata(M,N,delta,theta,f,SNR)
% Inputs --------------------------------------
% M: Number of antennas
% N: Number of samples per source
% delta: Antenna aperture in terms of wavelength
% theta: Directions of the sources in [-pi/2, pi/2)
% f: Normalized source frequencies in [0, 1)
% SNR: Desired noise to signal power ratio in dB
% Outputs --------------------------------------
% X: MxN received noisy sample matrix
% A: Mxd antenna response matrix
% S: dxN source samples

d = length(theta);          % Number of sources
k = 0:N-1;                  % Sample indices

% Building S
S = exp(1j* 2*pi*f(:)*k);

% Building A
A = zeros(M, d);
for i = 1:d
    angle_rad = deg2rad(theta(i));
    tau = (0:M-1).'*delta*sin(angle_rad);  
    A(:, i) = exp(1j*2*pi*tau);
end

X_clean = A * S;
% Compute corresponding noise power 
signal_power = mean(abs(X_clean(:)).^2);
SNR = 10^(SNR / 10);
noise_power = signal_power / SNR;
% Generate complex gaussian with the noise power
noise = sqrt(noise_power/2) * (randn(M, N) + 1j * randn(M, N));
% Add noise
X = X_clean + noise;
end

