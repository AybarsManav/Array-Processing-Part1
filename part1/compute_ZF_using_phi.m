function W_theta_ZF = compute_ZF_using_phi(phi_est, N, X)
d = length(phi_est);
S_structured = zeros(d, N);
for k = 1:d
    for n = 1:N
        S_structured(k, n) = exp(1j * 2 * pi * phi_est(k) * (n - 1));
    end
end
% From S obtain beamformer to isolate A
W_theta_ZF = (S_structured * pinv(X))';
end

