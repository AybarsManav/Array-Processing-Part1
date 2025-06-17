clc
clear
close all


delta = 1/2;
theta_range = linspace(-90,90,1000);
f = [0.1 0.3];
theta = [-20,30];
SNR = 20;
M = 5;
N = 20;
antenna = linspace(1,M,M);

%% data generation

[X,A,S] = gendata(M,N,delta,theta,f,SNR);

[U, sigma, V] = svd(X);

sing = diag(sigma);
figure()
scatter(antenna,sing)

%% ESPRIT

d = length(theta);
theta_est = esprit(X,d);
%theta in degrees
theta_est = asind(theta_est / (2 * pi * delta));

%% data generation - Channel estimation
c = 1/sqrt(2);
N = 500;
s = zeros(1,N);
for i=1:N
    if randi([0,1])==1
        a = 1;
    else
        a = -1;
    end
    if randi([0,1])==1
        b = 1;
    else
        b = -1;
    end

    s(i) = a*c + 1i*b*c;
end
%s = [a+a*1i, a-a*1i, -a+a*1i, -a-a*1i, a+a*1i];
N = length(s);
sigma = 0;
P = 8;

x = gendata_conv(s,P,N,sigma);

%% creating X pairs

xr = reshape(x, P, N).';

first = reshape(xr(1:499, :).', P,N-1);
second = reshape(xr(2:500,:).',P,N-1);
X = [first;second];

%% Zero Forcing and Wiener

c = 1/sqrt(2);
N = 500;
s = zeros(1,N);
for i = 1:N
    a = 2 * randi([0,1]) - 1;  
    b = 2 * randi([0,1]) - 1;
    s(i) = a*c + 1i*b*c;
end

sigma = 0.5;
P = 4;

x = gendata_conv(s,P,N,sigma);
xr = reshape(x, P, N).';

first = reshape(xr(1:499, :).', P,N-1);
second = reshape(xr(2:500,:).',P,N-1);
X = [first;second];


% Zero forcing
A = [ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4), zeros(1,P); 
    zeros(1,P), ones(1,P/4), -ones(1,P/4), ones(1,P/4), -ones(1,P/4)].';
W1 = pinv(A);
S1 = pinv(A)*X;
s1 = zeros(1,N);
s1(1:499) = S1(1,:);
s1(500) = S1(2,end);


%Wiener

W2 = inv(A*A' + sigma* eye(2*P,2*P))*A;
S2 = W2'*X;
s2 = zeros(1,N);
s2(1:499) = S2(1,:);
s2(500) = S2(2,end);


%choose different delays
                                %s1_d0 = s1(1:P:end);
%s2_d0 = s2(1:p:end);


figure()
scatter(real(s1), imag(s1), 'filled');
hold on
scatter(real(s2), imag(s2), 'filled');

title('QPSK Constellation');
xlabel('In-phase (I)');
ylabel('Quadrature (Q)');
legend('Zero-forcing','Wiener Equalizer')
axis equal;
grid on;

all_outputs_zf = zeros(2*P, N-1);
all_outputs_wiener = zeros(2*P, N-1);

for delay = 1:2*P

    
    
    w_zf_delay = W1(delay,:)
    w_wiener_delay =

    all_outputs_zf(delay, :) = w_zf_delay * X;        
    all_outputs_wiener(delay, :) = w_wiener_delay * X;
end

%% Plot QPSK constellation for each delay
figure;
for delay = 1:2*P
    subplot(4,4,delay)
    scatter(real(all_outputs_zf(delay, :)), imag(all_outputs_zf(delay, :)), 10, 'b', 'filled');
    hold on
    scatter(real(all_outputs_wiener(delay, :)), imag(all_outputs_wiener(delay, :)), 10, 'r', 'filled');
    title(['Delay ', num2str(delay)]);
    axis equal;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    grid on;
end
sgtitle('Constellations for Different Delays (Blue: ZF, Red: Wiener)');