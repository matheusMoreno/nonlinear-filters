function x_dot = mackeyglass_eq(x_t, x_t_minus_tau, a, b, n)
% MACKEYGLASS_EQ returns dx/dt of Mackey-Glass delayed differential equation
    x_dot = -b*x_t + a*x_t_minus_tau/(1 + x_t_minus_tau^n);
end