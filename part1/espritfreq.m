function f = esprit_freq(X,d)

N = size(X, 2);
Z = [X(:, 1:end-1) X(:, 2:end)];        % Smoothing in the temporal domain
[U_z, Sigma_z, V_z] = svd(Z, "econ");   
V_H_z = V_z';

U_z = U_z(:, 1:d); 
Sigma_z = Sigma_z(1:d, 1:d);
V_H_z = V_H_z(1:d, :);
V_H_x = V_H_z(1:d, 1:N-1);
V_H_y = V_H_z(1:d, N:end);
Matrix = V_H_y * (V_H_x' * inv(V_H_x * V_H_x'));
[T, PSI] = eig(Matrix);
f = angle(diag(PSI./ abs(PSI))) / (2 * pi);
end

