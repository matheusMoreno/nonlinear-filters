a        = 0.2;     % value for a in eq (1)
b        = 0.1;     % value for b in eq (1)
n        = 10;
tau      = 20;		% delay constant in eq (1)
x0       = 1.2;		% initial condition: x(t=0)=x0
deltat   = 0.1;	    % time step size (which coincides with the integration step)
len      = 5000;	% total no. of samples, excluding the given initial condition

[X, T] = mackeyglass(a, b, n, tau, x0, deltat, len);

figure
plot(T, X);
set(gca,'xlim',[0, T(end)]);
xlabel('t');
ylabel('x(t)');
title(sprintf('A Mackey-Glass time serie (tau=%d)', tau));