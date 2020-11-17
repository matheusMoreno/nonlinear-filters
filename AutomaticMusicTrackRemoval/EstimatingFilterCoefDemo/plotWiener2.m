%% Plotting results in dB
figure;

H = fft(filterCoeff,nfft);
plot(w(1:nfft/2+1)/pi*(Fs/2000),20*log10(abs(H(1:nfft/2+1))), 'b')
xlabel('Normalized Frequency (kHz)', 'Interpreter', 'LaTex');
ylabel('Magnitude (dB)', 'Interpreter', 'LaTex');
hold on 

%% WienerCoeff Plot
H = freqResp.averageCoeff;
plot(w(1:nfft/2+1)/pi*(Fs/2000),20*log10(abs(H(1:nfft/2+1))), 'r')

%% WienerCoeff2 Plot

H = freqResp.medianCoeff;
plot(w(1:nfft/2+1)/pi*(Fs/2000),20*log10(abs(H(1:nfft/2+1))), 'g')

H = freqResp.medianFreq;
plot(w(1:nfft/2+1)/pi*(Fs/2000),20*log10(abs(H(1:nfft/2+1))), 'c')

%% WienerCoeff2 Plot
%[H,w] = freqz(wienerCoeffMedian,1,512);

%subplot(2,1,1)
%H = wienerCoeffCluster_freq;
%plot(w(1:nfft/2+1)/pi,20*log10(abs(H(1:nfft/2+1))), '--c')
legend({'Ref', 'TDAV' , 'TDMV', 'FDMV'}, 'Interpreter', 'LaTex', 'Location','southwest');
%subplot(2,1,2)
%plot(w(1:nfft/2+1)/pi,unwrap(angle(H(1:nfft/2+1)))/pi, '--c')
%legend({'Ref', 'TDAV' , 'TDMV', 'FDMV', 'NK Clustering'}, 'Interpreter', 'LaTex');


%% Plotting results in linear scale
%figure;
%H = freqResp.medianFreqABS;
%[H,w] = freqz(amplitude*b,1,512); 
%plot(w/pi,(abs(H)))
% ax = gca;
% ax.YLim = [-100 10];
% %ax.XTick = 0:.25:1;
%xlabel('Normalized Frequency (\times \pi rad/sample)');
%ylabel('Magnitude');
%hold on
%[H,w] = freqz(wienerCoeff(:,1),1,512); 
%plot(w/pi,(abs(H)), 'r')
%[H,w] = freqz(wienerCoeff2,1,512); 
%plot(w/pi,(abs(H)), 'k')

