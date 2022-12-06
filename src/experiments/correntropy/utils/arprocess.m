function signal = arprocess(N, Ts, duration)
    %ARPROCESS generates a signal that behaves like a soundwave.

    assert(mod(N, 2) == 0, 'N has to be even.');

    N = N / 2;

    % Generate poles
    angles = pi * rand(N, 1);
    radius = 0.75 + (0.99 - 0.75) * rand(N, 1);
    poles = [radius .* exp(1i * angles); radius .* exp(-1i * angles)];

    % Create system
    sys = zpk([], poles, 1, Ts);

    % Generate desired signal (no normalization)
    t = 0:Ts:duration;
    u = wgn(length(t), 1, 0);
    signal = lsim(sys, u, t);

    % Normalize signal so that |s| <= 1
    maxVal = max(abs(signal));
    g = 1 / maxVal;
    signal = g .* signal;
    signal = signal(:);
end
