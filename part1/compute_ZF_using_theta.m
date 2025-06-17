function W_theta_ZF = compute_ZF_using_theta(theta_est, M)
% steering matrix from angle estimates
d = length(theta_est);
delta = 0.5;
A_theta = zeros(M, d);
for k = 1:d
    for m = 1:M
        A_theta(m,k) = exp(1j * 2 * pi * delta * sind(theta_est(k)) * (m - 1));
    end
end
W_theta_ZF = pinv(A_theta)'; 
end

