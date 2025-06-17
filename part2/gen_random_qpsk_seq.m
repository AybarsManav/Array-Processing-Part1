function s = gen_random_qpsk_seq(N)
% inputs:
% N: Number of QPSK symbols to generate
% outputs:
% s: qpsk data sequence  

% Generate N random bits for I (in‐phase) and Q (quadrature) branches
bitsI = randi([0 1], N, 1);    % each entry is 0 or 1
bitsQ = randi([0 1], N, 1);

% Map {0,1} to {−1,+1}
I = 2*bitsI - 1; 
Q = 2*bitsQ - 1;

% Form QPSK symbols and normalize to unit power:
s = (I + 1j*Q) / sqrt(2);
end