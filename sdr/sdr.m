function SDR = sdr(x, y)
%SDR Computes the signal-to-distortion ratio between x and y
%   This function returns the SDR between two signals x and y, considering
%   that y is a (possibly distorted) estimate of x. The SDR is defined in
%   Vincent et. al (2006) as
%
%   SDR = 10 * log10(||source||^2 / ||interf + noise + artif||^2)
%
%   Considering that y is an estimate with interference, noise and
%   atifacts, we can use x as a "true estimate" and everything else not
%   present in x but present in y can be considered the sum of those three
%   unwanted effects.

    SDR = 10 * log10(norm(x) ^ 2 / norm(x - y) ^ 2);
end