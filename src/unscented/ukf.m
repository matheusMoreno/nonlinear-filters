function d_hat = ukf(x, d, L, P, pNoise, mNoise)

    function w_next = stateFunW(w)
    % STATEFUN Function that updates the states of UKF model.
        w_next = w;
    end

    function y = measurementFunW(w, s, u)
    % MEASUREMENTFUNC Function that updates measurement of UKF model.
        coefs = [u; s; 1];
        for p = 1:P-1
            coefs = kron(coefs, [u; s; 1]);
        end
        y = w' * coefs;
    end

    function s_next = stateFunX(s, u)
    % STATEFUN Function that updates the states of UKF model.
        s_next = [u; s(1:end - 1)];
    end

    function y = measurementFunX(s, w, u)
    % MEASUREMENTFUNC Function that updates measurement of UKF model.
        y = measurementFunW(w, s, u);
    end


    assert(length(x) == length(d), 'The sizes of x and d must be the same.');

    x = x(:);
    d = d(:);
    N = length(d);

    % Initialize parameters
    w = zeros((L + 1) ^ P, 1);
    states_i = zeros(L - 1, 1);
    d_hat = zeros(N, 1);

    % Parameters filter
    ukf_w = unscentedKalmanFilter(@stateFunW, @measurementFunW);
    ukf_w.State = w;
    ukf_w.ProcessNoise = pNoise;
    ukf_w.MeasurementNoise = mNoise;

    % States filters
    ukf_x = unscentedKalmanFilter(@stateFunX, @measurementFunX);
    ukf_x.State = states_i;
    ukf_x.ProcessNoise = pNoise;
    ukf_x.MeasurementNoise = mNoise;

    for n = 1:N
        states = correct(ukf_x, d(n), w, x(n));
        w = correct(ukf_w, d(n), states, x(n));

        d_hat(n) = measurementFunX(states, w, x(n));

        predict(ukf_x, x(n));
        w = predict(ukf_w);

        if mod(n, 4000) == 0
            fprintf('At %d\n', n);
        end
    end
end
