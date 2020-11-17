clear all
close all;

[slong fs] = audioread('test.wav');
slong = 0.99*slong/max(abs(slong));

% Apply a distortion
alpha = 2;
xlong = (1/alpha)*atan(alpha*slong);
xlong = std(slong)*(xlong/std(xlong));
plot(slong,slong,'.');
hold on;
plot(slong,xlong,'.k');

% Apply the recovery algorithm
N = 1e3;
P = 5;
sigmaL0 = 1e-2;

stillong = duarteSmoothL0norm(xlong,N,P,sigmaL0);
stillong = std(slong)*(stillong/std(stillong));

% Scatter plots all signals
hold on
plot(slong,stillong,'.g');

