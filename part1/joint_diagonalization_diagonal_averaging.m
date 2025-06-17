function [f, theta] = joint_diagonalization_diagonal_averaging(X, m)
% joint_diagonalization_diagonal_averaging Performs time smoothing and joint diagonalization.
% Assumes d = 2 number of sources fixed!
%
%   [f, theta] = joint_diagonalization_diagonal_averaging(X, m)
%
%   Inputs:
%       X - Signal data matrix of size [M x N]
%       m - Smoothing factor 2 < m < M-1
%
%   Outputs:
%       f     - Estimated frequencies
%       theta - Estimated angular parameters (degrees)

d = 2; % Number of sources

[M, N] = size(X);

if m < 2 || m >= N - 1
    error('Parameter m must satisfy 1 < m < N - 1.');
end

% Time smoothing of spatially shifted data
X_m1 = time_smoothing(X(1:M-1, :), m);
X_m2 = time_smoothing(X(2:M, :), m);

% Concatenate for single spatial shift
X_cal = [X_m1; X_m2];

% SVD decomposition and truncation
[U, D, ~] = svd(X_cal, 'econ');

U = U(:, 1:d);
D = D(1:d, 1:d);

% Compute transition matrices
num_blocks = size(U, 1) / (M - 1);
M_matrices = zeros(d, num_blocks * d); % dxd square matrices of the form T^-1 D T

for i = 1:num_blocks/2-1
    rows = (i) * (M - 1) + 1 : (i+1) * (M - 1);
    pinv_U = pinv(U((i - 1) * (M - 1) + 1 : i * (M - 1),:));
    M_matrices(:,1:d) = eye(2,2);
    M_matrices(:, (i) * d + 1 : (i+1) * d) = pinv_U * U(rows, :);
end

for i = 1:num_blocks/2
    rows = (i-1) * (M - 1) +(2*(M-1)*m)/2 + 1 : (i) * (M - 1)+(2*(M-1)*m)/2;
    temp = (i - 1) * (M - 1) + 1 : i * (M - 1);
    pinv_U = pinv(U(temp,:));
    M_matrices(:, (i-1) * d + 1  + num_blocks: (i) * d + num_blocks) = pinv_U * U(rows, :);
end

% Perform joint diagonalization
[~, D_joint] = joint_diag(M_matrices, 1e-8);

% Frequency estimation
averaged_diagonal_phi = zeros(d, d);
for i=2:num_blocks/2
    averaged_diagonal_phi = averaged_diagonal_phi + D_joint(:, d * i - 1: d * i);
end
averaged_diagonal_phi = averaged_diagonal_phi / (num_blocks/2 - 1);
f = angle(diag(averaged_diagonal_phi) ./ abs(diag(averaged_diagonal_phi)) ) / (2 * pi);


% Angular parameter estimation
averaged_diagonal_theta = zeros(d, d);
lambda_start = num_blocks;
for i=1:num_blocks/2
    averaged_diagonal_theta = averaged_diagonal_theta + D_joint(:, lambda_start+ d * i - 1:lambda_start+ d * i);
end
averaged_diagonal_theta = averaged_diagonal_theta / (num_blocks/2);
lambda= diag(averaged_diagonal_theta ./ abs(averaged_diagonal_theta));
theta = asind(angle(lambda) / pi);

