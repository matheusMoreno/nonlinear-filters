fs = 2000;
Ts = 1 / fs;
duration = 1;
t = 0:Ts:duration;
t = t(:);
len = length(t);

% Sinewave with harmonicse
basef = 3;
As = [0.4 -0.2 -0.1 -0.2 0.4];
sine = zeros(length(t), 1);
for k = 1:length(As)
    sine = sine + As(k) .* sin(2 * pi * k * basef * t);
end
sine = sine(:);

sine_clipped = 2 .* atan(1.5 .* sine) ./ pi;
sine_clipped_n = sine_clipped + wgn(len, 1, -20);

% History size and number of polynomials
L = 3;
P = 4;

w = zeros((L + 1) ^ P, 1);
x_i = 5e-5 * rand(L - 1, 1);
sine_est = zeros(len, 1);

ukf_w = unscentedKalmanFilter(@stateFunW, @measurementFunW);
ukf_w.State = w;
ukf_w.ProcessNoise = 1e-2;
ukf_w.MeasurementNoise = 10;

ukf_x = unscentedKalmanFilter(@stateFunX, @measurementFunX);
ukf_x.State = x_i;
ukf_x.ProcessNoise = 1e-2;
ukf_x.MeasurementNoise = 10;


for k = 1:len
    u = sine(k);
    y = sine_clipped_n(k);

    x = correct(ukf_x, y, w, u);
    w = correct(ukf_w, y, x, u);

    sine_est(k) = measurementFunX(x, w, u);

    w = predict(ukf_w);
    x = predict(ukf_x, u);
end


disp(mean((sine_clipped - sine_est) .^ 2));

plot(t, sine_clipped)
hold on
%plot(t, sine)
plot(t, sine_est)
saveas(gcf, 'sine.png')
hold off


function w_next = stateFunW(w)
% STATEFUN Function that updates the states of UKF model.
    w_next = w;
end


function y = measurementFunW(w, x, u)
% MEASUREMENTFUNC Function that updates measurement of UKF model.
    u_x = [u; x; 1];
    P = round(log(length(w)) / log(length(u_x)));

    coefs = u_x;
    for p = 1:P-1
        coefs = kron(coefs, u_x);
    end

    y = w' * coefs;
end


function x_next = stateFunX(x, u)
% STATEFUN Function that updates the states of UKF model.
    x_next = [u; x(1:end - 1)];
end


function y = measurementFunX(x, w, u)
% MEASUREMENTFUNC Function that updates measurement of UKF model.
    u_x = [u; x; 1];
    P = round(log(length(w)) / log(length(u_x)));

    coefs = u_x;
    for p = 1:P-1
        coefs = kron(coefs, u_x);
    end

    y = w' * coefs;
end
