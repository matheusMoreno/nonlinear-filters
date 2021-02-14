% Test cases for the nonlinear discrete Wiener filter
% By Matheus F. Moreno
% SMT/DEL/UFRJ - Escola Politecnica

%% Instantiating the input signal: a discrete sine wave
x = sin((0:1000) * pi / 100);
x = x(:);

figure
hold on
plot(1:length(x), x);
set(gca,'xlim',[0, length(x)]);
xlabel('n');
ylabel('x(n)');
title('Original signal');


%% First case: filter estimation
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
sig_len = length(y_filt);

% Filter specifications
sizeFilters = 10;
hop = 500;
sigma = 10;

% Creating the filters
[y_filt_est, MSE] = kernelWiener(x, y_filt, sizeFilters, hop, sigma);
fprintf('MSE for filter estimation: %d\n', MSE);

% Plotting the filtered signal and the estimation
figure
hold on
plot(1:sig_len, y_filt);
plot(1:sig_len, y_filt_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Kernel Wiener Filter for filter estimation');


%% Second case: one step prediction
y_pred = x(2:end);
x_pred = x(1:end - 1);
sig_len = length(y_pred);

% Filter specifications
sizeFilters = 10;
hop = 200;
sigma = 0.16;

% Estimating the signal
[y_pred_est, MSE] = kernelWiener(x_pred, y_pred, sizeFilters, hop, sigma);
fprintf('MSE for one step prediction: %d\n', MSE);

figure
hold on
plot(1:sig_len, y_pred);
plot(1:sig_len, y_pred_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Kernel Wiener Filter for one step prediction');


%% Third case: noisy filtered signal
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
y_noisy = y_filt + randn(length(y_filt), 1) * 0.1;
sig_len = length(y_filt);

% Filter specifications
sizeFilters = 10;
hop = 200;
sigma = 3;

% Estimating the signal
[y_noisy_est, ~] = kernelWiener(x, y_noisy, sizeFilters, hop, sigma);
MSE = mean((y_filt - y_noisy_est) .^ 2);
fprintf('MSE for filtered noisy signal estimation: %d\n', MSE);

figure
subplot(1,2,1);
hold on
plot(1:sig_len, x);
plot(1:sig_len, y_noisy);
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('x(n), y(n) + w(n)');
title('Filtered signal with noise');
subplot(1,2,2)
hold on
plot(1:sig_len, y_filt);
plot(1:sig_len, y_noisy_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Kernel Wiener Filter for filtered signal noise reduction');


%% Fourth case: nonlinear filter
y_nonl = zeros(length(x), 1);
sig_len = length(x);
x0 = 0;
x_m1 = 0;
y0 = 0;

% Creating a very weird nonlinear filter
for i = 1:sig_len,
    y_nonl(i) = sin(abs(x(i))^(1/2)) * exp(x0 + x_m1) + cos(y0);
    x_m1 = x0;
    x0 = x(i);
    y0 = y_nonl(i);
end;

% Filter specifications
sizeFilters = 10;
hop = 100;
sigma = 0.05;

% Estimating the signal
[y_nonl_est, MSE] = kernelWiener(x, y_nonl, sizeFilters, hop, sigma);
fprintf('MSE for nonlinear prediction: %d\n', MSE);

figure
hold on
plot(1:sig_len, y_nonl);
plot(1:sig_len, y_nonl_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Kernel Wiener Filter for nonlinear prediction');


%% Fifth case: noisy nonlinear filter
y_nnl = zeros(length(x), 1);
N = length(x);
x0 = 0;
x_m1 = 0;
y0 = 0;

% Creating a very weird nonlinear filter
for i = 1:N,
    y_nnl(i) = sin(abs(x(i))^(1/2)) * exp(x0 + x_m1) + cos(y0);
    x_m1 = x0;
    x0 = x(i);
    y0 = y_nnl(i);
end;

y_nonl = y_nnl;
y_nnl = y_nnl + 0.5 * randn(N, 1);

% Filter specifications
sizeFilters = 10;
hop = 200;
sigma = 0.05;

% Estimating the signal
[y_nonl_est, ~] = kernelWiener(x, y_nnl, sizeFilters, hop, sigma);
MSE = mean((y_nonl - y_nonl_est) .^ 2);
fprintf('MSE for nonlinear noisy prediction: %d\n', MSE);

figure
subplot(1,2,1);
hold on
plot(1:N, x);
plot(1:N, y_nnl);
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('x(n), y(n) + w(n)');
title('Filtered signal with noise');
subplot(1,2,2)
hold on
plot(1:N, y_nonl);
plot(1:N, y_nonl_est, '--');
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Nonlinear Wiener Filter for nonlinear noisy prediction');