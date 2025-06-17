function [f, theta] = joint_diagonalization(X, m)
% JOINT_DIAGONALIZATION Performs time smoothing and joint diagonalization.
% Assumes d = 2 number of sources fixed!
%
%   [f, theta] = JOINT_DIAGONALIZATION(X, m)
%
%   Inputs:
%       X - Signal data matrix of size [M x N]
%       m - Smoothing factor 2 < m < M-1
%
%   Outputs:
%       f     - Estimated frequencies
%       theta - Estimated angular parameters (degrees)

d = 2; % Number of sources
m = 15;

[M, N] = size(X);

if m <= 2 || m >= N - 1
    error('Parameter m must satisfy 2 < m < N - 1.');
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
 f_arr = zeros(floor((num_blocks)/2),2);
 for i=1:num_blocks/2
     f_arr(i,:) = angle(diag(D_joint(:, d * i - 1: d * i) ./ abs(D_joint(:, d * i - 1: d * i)))) / (2 * pi);
 end

f = mean(f_arr(2:end,:));
theta_arr = zeros((num_blocks)/2,2);
% Angular parameter estimation
lambda_start = num_blocks;
for i=1:num_blocks/2
    lambda= diag(D_joint(:, lambda_start+ d * i - 1:lambda_start+ d * i) ./ ...
              abs(D_joint(:, lambda_start+d * i - 1:lambda_start+ d * i)));
    theta_arr(i,:) = asind(angle(lambda) / pi);
end
theta_arr = theta_arr(2:end,:);
theta = mean(theta_arr);


