function theta = esprit(X,d)

[M, N] = size(X);
[U, sigma, V] = svd(X);

Us = U(:, 1:d);  

Ux = Us(1:M-1, :);
Uy = Us(2:M, :);

est = pinv(Ux)*Uy;
[V, D] = eig(est);
lambda = diag(D ./ abs(D));

theta = angle(lambda);
theta = asind(theta / pi);
end

