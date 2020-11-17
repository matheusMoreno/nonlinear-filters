function [X, T] = mackeyglass(a, b, n, tau, x0, dt, len)
%MACKEYGLASS returns a Mackey-Glass time series of type 
% x' = -b*x(t) + a*x(t-tau)/(1 + x(t-tau)^n), initial condition x(0) = x0,
% step dt and size (number of samples) 'len'
%
% Based on this code:
% https://www.mathworks.com/matlabcentral/fileexchange/24390-mackey-glass-time-series-generator
%

    time = 0;
    ind_t_minus_tau = 1;
    init_length = floor(tau/dt); % x(t) = 0 for all t < 0
    x_t = x0;

    X = [zeros(init_length, 1); zeros(len + 1, 1)]; % vector of samples
    T = zeros(len + 1, 1); % vector of time samples

    for i = 1:len + 1,
        % Update values for current iteration
        X(init_length + i) = x_t;
        T(i) = time;
        x_t_minus_tau = X(ind_t_minus_tau);
  
        % Calculate the current value
        x_t_plus_deltat = mackeyglass_rk4(x_t, x_t_minus_tau, dt, a, b, n);

        % Update values for next iteration
        ind_t_minus_tau = ind_t_minus_tau + 1;
        time = time + dt;
        x_t = x_t_plus_deltat;
    end

    % Remove the initial values
    X = X(init_length + 1:end);
end