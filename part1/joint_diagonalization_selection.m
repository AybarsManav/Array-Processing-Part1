function [f, theta] = joint_diagonalization_selection(X, m)
% joint_diagonalization_selection Performs time smoothing and joint diagonalization.
% Assumes d = 2 number of sources fixed!
%
%   [f, theta] = joint_diagonalization_selection(X, m)
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

% Frequency
U_x_phi = U(1:(M-1)*m-(M-1), :);        % First (m-1) * (M-1) rows
U_y_phi = U(M:(M-1)*m, :);              % Last (m-1) * (M-1) rows
M_phi = pinv(U_x_phi) * U_y_phi;

% Theta
U_x_theta = U(1:(M-1)*m, :);                % First half
U_y_theta = U((M-1)*m + 1: 2*(M-1)*m, :);   % Last half
M_theta = pinv(U_x_theta) * U_y_theta;

M_matrices = [M_phi, M_theta];

% Perform joint diagonalization
[~, D_joint] = joint_diag(M_matrices, 1e-8);

D_freq_diag = diag(D_joint(:, 1:d));
D_theta_diag = diag(D_joint(:, d+1:2*d));

f = angle(D_freq_diag ./ abs(D_freq_diag)) / (2 * pi);
theta = asind(angle(D_theta_diag ./ abs(D_theta_diag)) / pi );

end


