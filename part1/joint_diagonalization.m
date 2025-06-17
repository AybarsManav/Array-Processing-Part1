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
    U_11 = U(1:(M-1), :);
    pinv_U_11 = pinv(U_11);
    num_blocks = size(U, 1) / (M - 1);
    M_matrices = zeros(d, num_blocks * d); % dxd square matrices of the form T^-1 D T

    for i = 1:num_blocks
        rows = (i - 1) * (M - 1) + 1 : i * (M - 1);
        M_matrices(:, (i - 1) * d + 1 : i * d) = pinv_U_11 * U(rows, :);
    end

    % Perform joint diagonalization
    [~, D_joint] = joint_diag(M_matrices, 1e-8);

    % Frequency estimation
    f = angle(diag(D_joint(:, 3:4) ./ abs(D_joint(:, 3:4)))) / (2 * pi);

    % Angular parameter estimation
    lambda_start = num_blocks;
    lambda = diag(D_joint(:, lambda_start + 1 : lambda_start + 2) ./ ...
                  abs(D_joint(:, lambda_start + 1 : lambda_start + 2)));
    theta = asind(angle(lambda) / pi);
end
