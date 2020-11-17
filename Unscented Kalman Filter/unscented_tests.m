% Test cases for the Unscented Kalman Filter
% By Matheus F. Moreno - 10/2020
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


%% UKF hyperparameters
%  They'll stay the same throughout the entire testing, for simplicity

n = 6;                                      % Delay size (max = n - 1)
p = 4;                                      % Number of terms in the poly

q = 0.1;                                    % "Inovation" for parameters
r = 0.5;                                    % Std of measurement
Q = q^2 * eye(n * p);                       % Covariance of process
R = r^2;                                    % Covariance of measurement

f = @(x)x;                                  % State equations


%% First case: filter estimation
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
N = length(y_filt);

% Creating the polynomial coefficients from the original signal
u_buffer = buffer(x, n, n - 1);
coefs = zeros(n * p, N);
for j = 0:p - 1,
    coefs(j*n + 1:j*n + n, :) = u_buffer .^ (2 * j + 1);
end;

s = repmat([zeros(n - 1, 1); 1], p, 1);     % Initial state
u = s + q * randn;                          % Initial state with noise
P = eye(n * p);                             % Initial state covraiance
uV = zeros(n * p, N);                       % Estimate for parameters
sV = zeros(n * p, N);                       % Actual parameters

% Estimating y_filt
y_filt_est = zeros(N, 1);
for k = 1:N,
    h = @(x)x' * coefs(:, k);               % Defining the observation function
    z = y_filt(k);                          % Measurement
    sV(:, k) = s;                           % Actual state (parameters)
    [u, P] = ukf(f, u, P, h, z, Q, R);      % ukf
    uV(:, k) = u;                           % Save estimate
    s = f(s) + q * randn(n * p, 1);         % Update process
    y_filt_est(k) = h(u);                   % Estimate filter value
end;

% Creating the filters
MSE = mean((y_filt - y_filt_est) .^ 2);
fprintf('MSE for filter estimation: %d\n', MSE);

% Plotting the filtered signal and the estimation
figure
hold on
plot(1:N, y_filt);
plot(1:N, y_filt_est, '--');
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Unscented Kalman Filter for filter estimation');


%% Second case: one step prediction
y_pred = x(2:end);
x_pred = x(1:end - 1);
N = length(y_pred);

% Creating the polynomial coefficients from the original signal
u_buffer = buffer(x_pred, n, n - 1);
coefs = zeros(n * p, N);
for j = 0:p - 1,
    coefs(j*n + 1:j*n + n, :) = u_buffer .^ (2 * j + 1);
end;

s = repmat([zeros(n - 1, 1); 1], p, 1);     % Initial state
u = s + q * randn;                          % Initial state with noise
P = eye(n * p);                             % Initial state covraiance
uV = zeros(n * p, N);                       % Estimate for parameters
sV = zeros(n * p, N);                       % Actual parameters

% Estimating y_pred
y_pred_est = zeros(N, 1);
for k = 1:N,
    h = @(x)x' * coefs(:, k);               % Defining the observation function
    z = y_pred(k);                          % Measurement
    sV(:, k) = s;                           % Actual state (parameters)
    [u, P] = ukf(f, u, P, h, z, Q, R);      % ukf
    uV(:, k) = u;                           % Save estimate
    s = f(s) + q * randn;                   % Update process
    y_pred_est(k) = h(u);                   % Estimate filter value
end;

% Creating the filters
MSE = mean((y_pred - y_pred_est) .^ 2);
fprintf('MSE for one step prediction: %d\n', MSE);

figure
hold on
plot(1:N, y_pred);
plot(1:N, y_pred_est, '--');
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Unscented Kalman Filter for one step prediction');


%% Third case: noisy filtered signal
filt = [-3, -2, -1, 1, 2, 3];
y_filt = fftfilt(filt, x);
N = length(y_filt);
y_noisy = y_filt + r * randn(N, 1);

% Creating the polynomial coefficients from the original signal
u_buffer = buffer(x, n, n - 1);
coefs = zeros(n * p, N);
for j = 0:p - 1,
    coefs(j*n + 1:j*n + n, :) = u_buffer .^ (2 * j + 1);
end;

s = repmat([zeros(n - 1, 1); 1], p, 1);     % Initial state
u = s + q * randn;                          % Initial state with noise
P = eye(n * p);                             % Initial state covraiance
uV = zeros(n * p, N);                       % Estimate for parameters
sV = zeros(n * p, N);                       % Actual parameters

% Estimating y_noisy
y_noisy_est = zeros(N, 1);
for k = 1:N,
    h = @(x)x' * coefs(:, k);               % Defining the observation function
    z = y_noisy(k);                         % Measurement
    sV(:, k) = s;                           % Actual state (parameters)
    [u, P] = ukf(f, u, P, h, z, Q, R);      % ukf
    uV(:, k) = u;                           % Save estimate
    s = f(s) + q * randn;                   % Update process
    y_noisy_est(k) = h(u);                  % Estimate filter value
end;

MSE = mean((y_filt - y_noisy_est) .^ 2);
fprintf('MSE for noisy filtered signal: %d\n', MSE);

figure
subplot(1,2,1);
hold on
plot(1:N, x);
plot(1:N, y_noisy);
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('x(n), y(n) + w(n)');
title('Filtered signal with noise');
subplot(1,2,2)
hold on
plot(1:N, y_filt);
plot(1:N, y_noisy_est, '--');
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Unscented Kalman Filter for filtered signal noise reduction');


%% Fourth case: nonlinear filter
y_nonl = zeros(length(x), 1);
N = length(x);
x0 = 0;
x_m1 = 0;
y0 = 0;

% Creating a very weird nonlinear filter
for i = 1:N,
    y_nonl(i) = sin(abs(x(i))^(1/2)) * exp(x0 + x_m1) + cos(y0);
    x_m1 = x0;
    x0 = x(i);
    y0 = y_nonl(i);
end;

% Creating the polynomial coefficients from the original signal
u_buffer = buffer(x, n, n - 1);
coefs = zeros(n * p, N);
for j = 0:p - 1,
    coefs(j*n + 1:j*n + n, :) = u_buffer .^ (2 * j + 1);
end;

s = repmat([zeros(n - 1, 1); 1], p, 1);     % Initial state
u = s + q * randn;                          % Initial state with noise
P = eye(n * p);                             % Initial state covraiance
uV = zeros(n * p, N);                       % Estimate for parameters
sV = zeros(n * p, N);                       % Actual parameters

% Estimating y_nonl
y_nonl_est = zeros(N, 1);
for k = 1:N,
    h = @(x)x' * coefs(:, k);               % Defining the observation function
    z = y_nonl(k);                          % Measurement
    sV(:, k) = s;                           % Actual state (parameters)
    [u, P] = ukf(f, u, P, h, z, Q, R);      % ukf
    uV(:, k) = u;                           % Save estimate
    s = f(s) + q * randn;                   % Update process
    y_nonl_est(k) = h(u);                   % Estimate filter value
end;

MSE = mean((y_nonl - y_nonl_est) .^ 2);
fprintf('MSE for nonlinear prediction: %d\n', MSE);

figure
hold on
plot(1:N, y_nonl);
plot(1:N, y_nonl_est, '--');
set(gca,'xlim',[0, N]);
xlabel('n');
ylabel('y(n), y_est(n)');
title('Unscented Kalman Filter for nonlinear prediction');


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
y_nnl = y_nnl + r * randn(N, 1);

% Creating the polynomial coefficients from the original signal
u_buffer = buffer(x, n, n - 1);
coefs = zeros(n * p, N);
for j = 0:p - 1,
    coefs(j*n + 1:j*n + n, :) = u_buffer .^ (2 * j + 1);
end;

s = repmat([zeros(n - 1, 1); 1], p, 1);     % Initial state
u = s + q * randn;                          % Initial state with noise
P = eye(n * p);                             % Initial state covraiance
uV = zeros(n * p, N);                       % Estimate for parameters
sV = zeros(n * p, N);                       % Actual parameters

% Estimating y_nonl
y_nonl_est = zeros(N, 1);
for k = 1:N,
    h = @(x)x' * coefs(:, k);               % Defining the observation function
    z = y_nnl(k);                           % Measurement
    sV(:, k) = s;                           % Actual state (parameters)
    [u, P] = ukf(f, u, P, h, z, Q, R);      % ukf
    uV(:, k) = u;                           % Save estimate
    s = f(s) + q * randn;                   % Update process
    y_nonl_est(k) = h(u);                   % Estimate filter value
end;

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
title('Unscented Kalman Filter for nonlinear noisy prediction');
