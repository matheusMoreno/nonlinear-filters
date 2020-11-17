function x_t_plus_deltat = mackeyglass_rk4(x_t, x_t_minus_tau, deltat, a, b, n)
% MACKEYGLASS_RK4 computes the numerical solution of the Mackey-Glass delayed
% differential equation using the 4-th order Runge-Kutta method
    k1 = deltat*mackeyglass_eq(x_t,          x_t_minus_tau, a, b, n);
    k2 = deltat*mackeyglass_eq(x_t+0.5*k1,   x_t_minus_tau, a, b, n);
    k3 = deltat*mackeyglass_eq(x_t+0.5*k2,   x_t_minus_tau, a, b, n);
    k4 = deltat*mackeyglass_eq(x_t+k3,       x_t_minus_tau, a, b, n);
    x_t_plus_deltat = (x_t + k1/6 + k2/3 + k3/3 + k4/6);
end