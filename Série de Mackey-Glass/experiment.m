%% Creating the Mackey-Glass series
% Setting the parameters:
% x' = -b*x(t) + a*x(t-tau)/(1 + x(t-tau)^n), x0, step deltat, size len
a = 0.2;
b = 0.1;
n = 10;
tau = 20;
x0 = 1.2;
deltat = 0.1;
len = 2000;

% Generating the Mackey-Glass series
[S, T] = mackeyglass(a, b, n, tau, x0, deltat, len);

%% Testing the linear filter
% Setting the filter parameters
sizeFilters = 10;
hop = 100;

% Delayed input and output
input = S(1:end - 1);
output = S(2:end);

% Creating the filters
[y_est, ~, MSE] = wiener(input, output, sizeFilters, hop);
disp(MSE);

%x_est = nonlinearWiener(input, output, sizeFilters, 0.1);

% Estimating the values
%numEstimations = size(estimatedFilters, 2);
%xBuffer = buffer(input, sizeFilters, hop);
%x_est = [x0; zeros(numEstimations, 1)];

%for i = 1:numEstimations,
%    x_est(i + 1) = estimatedFilters(i)' * flip(xBuffer(i));
%end;

%% Plotting the results
figure
hold on
plot(T, S);
plot(T(2:end), y_est);
set(gca,'xlim',[0, T(end)]);
xlabel('t');
ylabel('x(t)');
title(sprintf('A Mackey-Glass time serie (tau=%d)', tau));