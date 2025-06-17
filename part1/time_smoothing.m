function [X_m] = time_smoothing(X, m)
M = size(X, 1);
N = size(X, 2);
window_size = N - m + 1;
X_m = zeros(m * M, window_size);
for i = 1:m
    X_m(1 + (i - 1)* M: i * M, :) = X(:, i:window_size + (i - 1));
end
end

