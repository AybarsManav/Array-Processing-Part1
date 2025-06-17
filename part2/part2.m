% Generate a random QPSK data symbol sequence
N = 500;                                    % Number of QPSK data symbols
sigma = 0.5;                                % Std of zero mean complex gaussian noise
P = 4;                                      % Upsampling factor
s = gen_random_qpsk_seq(N);                 % QPSK data symbols

x = gendata_conv(s, P, N, sigma);           % Received samples
X = reshape(x, P, []);                      % Intermediate data matrix where each column comes from a symbol
cal_X = [X(:, 1:end-1); X(:, 2:end)];       % Data matrix defined in (1) 2P x (N - 1)

% We can obtain cal_X as outer product of 2 rank-2 matrices, H (channel conv) and S
% cal_X = H * S
H = [ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4), zeros(1,P); 
    zeros(1,P), ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4)].';

% Zero forcing beamformer
W_ZF_H = pinv(H);                               % Zero-forcing beamformer
S_hat_ZF = W_ZF_H*cal_X;
s_hat_zf = zeros(1,N);
s_hat_zf(1:499) = S_hat_ZF(1,:);
s_hat_zf(500) = S_hat_ZF(2,end);

% Wiener beamformer
Rx_hat = cal_X * cal_X' / N;                    % Sample data covariance matrix
% Since E[SS'] = N * I
Rxs_hat = H;                                    % Sample correlation between symbols and received data
W_wiener = inv(Rx_hat) * Rxs_hat;               % If there is noise R_x is invertible
% W_wiener = inv(H*H' + sigma* eye(2*P,2*P))*H;
S_hat_wiener = W_wiener'*cal_X;                 % Reconstruct S using wiener
s_hat_wiener = zeros(1,N);
s_hat_wiener(1:499) = S_hat_wiener(1,:);
s_hat_wiener(500) = S_hat_wiener(2,end);

% Choose a delay (which row of W^H to take)
chosen_row = 1;
w_zf_H = W_ZF_H(chosen_row, :);
w_wiener_H = W_wiener(:, chosen_row)';

estimate_zf = w_zf_H * cal_X;
estimate_wiener = w_wiener_H * cal_X;

% Plotting
figure;
scatter(real(estimate_zf), imag(estimate_zf), 10, 'b', 'filled');
hold on
scatter(real(estimate_wiener), imag(estimate_wiener), 10, 'r', 'filled');
axis equal;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
xlabel("Real");
ylabel("Imag")
grid on;

sgtitle('Constellations for beamformers (Blue: ZF, Red: Wiener)');


%% Second part with doubled P
% Generate a random QPSK data symbol sequence
N = 500;                                    % Number of QPSK data symbols
sigma = 0.5;                                % Std of zero mean complex gaussian noise
P = 8;                                      % Upsampling factor
s = gen_random_qpsk_seq(N);                 % QPSK data symbols

x = gendata_conv(s, P, N, sigma);           % Received samples
X = reshape(x, P, []);                      % Intermediate data matrix where each column comes from a symbol
cal_X = [X(:, 1:end-1); X(:, 2:end)];       % Data matrix defined in (1) 2P x (N - 1)

% We can obtain cal_X as outer product of 2 rank-2 matrices, H (channel conv) and S
% cal_X = H * S
H = [ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4), zeros(1,P); 
    zeros(1,P), ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4)].';

% Zero forcing beamformer
W_ZF_H = pinv(H);                               % Zero-forcing beamformer
S_hat_ZF = W_ZF_H*cal_X;
s_hat_zf = zeros(1,N);
s_hat_zf(1:499) = S_hat_ZF(1,:);
s_hat_zf(500) = S_hat_ZF(2,end);

% Wiener beamformer
Rx_hat = cal_X * cal_X' / N;                    % Sample data covariance matrix
% Since E[SS'] = N * I
Rxs_hat = H;                                    % Sample correlation between symbols and received data
W_wiener = inv(Rx_hat) * Rxs_hat;               % If there is noise R_x is invertible
% W_wiener = inv(H*H' + sigma* eye(2*P,2*P))*H;
S_hat_wiener = W_wiener'*cal_X;
s_hat_wiener = zeros(1,N);
s_hat_wiener(1:499) = S_hat_wiener(1,:);
s_hat_wiener(500) = S_hat_wiener(2,end);

% Choose a delay (which row of W^H to take)
chosen_row = 1;
w_zf_H = W_ZF_H(chosen_row, :);
w_wiener_H = W_wiener(:, chosen_row)';

estimate_zf = w_zf_H * cal_X;
estimate_wiener = w_wiener_H * cal_X;

figure;
scatter(real(estimate_zf), imag(estimate_zf), 10, 'b', 'filled');
hold on
scatter(real(estimate_wiener), imag(estimate_wiener), 10, 'r', 'filled');
axis equal;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
xlabel("Real");
ylabel("Imag")
grid on;

sgtitle('Constellations for beamformers (Blue: ZF, Red: Wiener)');