% Test cases for the linear discrete Wiener filter
% By Matheus F. Moreno
% SMT/DEL/UFRJ - Escola Politecnica

%% Instantiating the input signal: a discrete sine wave
x = sin((0:2000) * pi / 100);
x = x(:);

figure
hold on
plot(1:length(x), x);
set(gca,'xlim',[0, length(x)]);
xlabel('n');
ylabel('x(n)');
title('Original signal');


%% First case: filter estimation with wiener()
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
sig_len = length(y_filt);

% Filter specifications
size_filter = 10;

% Creating the filters
y_filt_est = wiener(x, y_filt, size_filter);

% Plotting the filtered signal and the estimation
figure
hold on
plot(1:sig_len, y_filt);
plot(1:sig_len, y_filt_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Linear Wiener filter for filter estimation');


%% Second case: one step prediction with wienerCola()
y_pred = x(2:end);
x_pred = x(1:end - 1);
sig_len = length(y_pred);

% Filter specifications
size_filter = 20;
hop = 100;

% Estimating the signal
y_pred_est = wienerCola(x_pred, y_pred, size_filter, hop);

figure
hold on
plot(1:sig_len, y_pred);
plot(1:sig_len, y_pred_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Linear Wiener filter for one step prediction');


%% Third case: noisy filtered signal with wienerCola()
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
y_noisy = y_filt + randn(length(y_filt), 1) * 0.25;
sig_len = length(y_filt);

% Filter specifications
size_filter = 34;
hop = 256;

% Estimating the signal
y_noisy_est = wienerCola(x, y_noisy, size_filter, hop);

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
title('Linear Wiener filter for filtered signal noise reduction');


%% Fourth case: nonlinear filter with wienerCola()
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
size_filter = 10;
hop = 50;

% Estimating the signal
y_nonl_est = wienerCola(x, y_nonl, size_filter, hop);

figure
hold on
plot(1:sig_len, y_nonl);
plot(1:sig_len, y_nonl_est, '--');
set(gca,'xlim',[0, sig_len]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Linear Wiener filter for nonlinear prediction');

%% Fifth case: noisy nonlinear filter with wienerCola()
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
size_filter = 10;
hop = 50;

% Estimating the signal
y_nonl_est = wienerCola(x, y_nnl, size_filter, hop);

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
title('Linear Wiener Filter for nonlinear noisy prediction');