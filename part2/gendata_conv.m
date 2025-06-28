function x = gendata_conv(s,P,N,sigma)
% Generate received data matrix given:
% s: symbol vector
% P: Upsampling rate
% sigma: standard deviation of complex zero mean gaussian noise
% ------------------
% returns:
% X: Received samples matrix

    h = [ones(1, P/4), -ones(1, P/4), ones(1, P/4), -ones(1, P/4)];
    s_up = zeros(1, N*P);
    s_up(1:P:end) = s;
    
    x_clean = conv(s_up, h, "full");
    x_clean = x_clean(1:end - (P - 1));
    
    noise = sigma.^2/sqrt(2) * (randn(size(x_clean)) + 1j*randn(size(x_clean)));
    
    x = x_clean + noise;
    
end

