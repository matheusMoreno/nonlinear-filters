%close all;
clear all; 

%% Getting some important parameters for soundtrack removal
SetParametersEstGain;

 %% Reading Dialogue and Soundtrack Files
[dialogue, ~] = audioread(DATA_FILE);
dialogue = mean(dialogue, 2);
[soundtrack, Fs] = audioread(SOUNDTRACK_FILE);
soundtrack = mean(soundtrack, 2);
 
excerpt = soundtrack(initialSample:finalSample);   % the excerpt of the soundtrack that is going to be added to the dialogue signal

%% Applying gain on the excerpt before adding it to the dialogue signal

if variableGain_bool
    [gainCurve] = GenerateGainCurve(excerpt, position_vector, variableGain_vector);
else
    [gainCurve] = GenerateGainCurve(excerpt, constantGainPosition_vector, constantGain_vector);
end
 excerptEqualized = excerpt.*gainCurve;
 
%% Use this for a sinusoidal gain Curve
T = length(gainCurve);
%T1 = round(T/3);
%T2 = round(2*T/3);

gainSine = chirp([0:1/Fs:(T-1)/Fs], 0, (T-1)/Fs , 0.7)/2 + 0.7;
excerptEqualized = excerpt.*gainSine';   % comment this to not use the sinusoidal gain curve.

%% Use this for exponeltial + silent + constant gain curve
T1 = round(5*Fs);
T2 = round(8*Fs);
T3 = round(11*Fs);

gainExp(1:T1) = 0.2.^([1:T1]./Fs)+ 0.4;  % exponential part until 0.4 gain
gainExp(T1+1 : T2) = 0;                  % silence
gainExp(T2+1 : T3) = 0.8;                % constant = 0.8
gainExp(T3+1 : T ) = 0.5;                % constant = 0.5
%excerptEqualized = excerpt.*gainExp';    % comment this to not use an exponencial gain.

%% Creating the mixture file
[mixture, mix] = CreateMixture(dialogue, excerptEqualized, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue
 

%% Estimating the Gain Curve used in the Mixture Process

estimatedDelay = delay;   % supposing we have a way to estimate the delay correctly

parameters = {};
parameters.gainWindow = gainWindow;
parameters.gainHop = gainHop;
    
estGainCurve = EstimateGain(mix, excerpt, parameters);
 
 %mixture(delay:1:delay+length(excerpt) - 1) = mix - excerpt_low.*estGainCurve;
 
%  nloops = ceil((length(excerpt) - gainWindow) / gainHop) + 1;
%  
%  alpha = zeros(1,nloops);
%  alpha_0 = zeros(1,nloops);
%  mean_t = zeros(1,nloops);
%  
%  %alpha_time = estimatedDelay + floor(sizeWindow/2) + 1 +[0:nloops-1].*hop;
%  alpha_time = estimatedDelay + [0:nloops-1].*gainHop;
%  
%  i = 0;
%  while i < nloops
%     if i == nloops - 1
%         numberZeros = max(0,i*gainHop + gainWindow - length(excerpt));    
%         excerpt = [excerpt; zeros(numberZeros,1)];
%         numberZeros = max(0,estimatedDelay + i*gainHop + gainWindow - 1 - length(mixture));
%         mixture = [mixture; zeros(numberZeros,1)];
%         
%     end
%     t = excerpt( 1 + i*gainHop : i*gainHop + gainWindow); 
%     s = mixture(estimatedDelay + i*gainHop : estimatedDelay + i*gainHop + gainWindow - 1);
%     
%     s_0 = s - mean(s);       % zero mean
%     t_0 = t - mean(t);       % zero mean
%     mean_t(i+1) = mean(t);
%     alpha(i+1) = (s'*t)/(t'*t);          % estimated amplitude
%     alpha_0(i+1) = (s_0'*t_0)/(t_0'*t_0);  % estimated amplitude zero mean
%     if i > 0
%         
%     end    
%     i = i+1;
%  end
% 
%  nloops = length(alpha_time);
% for i = 1:nloops
%    if i == nloops
%        gainCurve(position_vector(i) : end) = gain_vector(i);
%    else
%        gainCurve(position_vector(i) : position_vector(i+1) - 1) = linspace(gain_vector(i), gain_vector(i+1), position_vector(i+1) - position_vector(i));        
%    end 
%  

%figure;
%plot(([0:length(gainCurve) - 1] + delay )/Fs, gainCurve);
%hold on
%plot(([10000:length(estGainCurve) - 1] + estimatedDelay )/Fs, estGainCurve(10001:end));

%plot(alpha_time/Fs, alpha, 'k');
%plot(alpha_time/Fs, alpha_0, 'r');
%title('Soundtrack Gain Estimation')
%xlabel('Time (s)');
%ylabel('Amplitude Gain');
%legend('Original' , 'Estimated');
%legend('Original' , 'Estimated', 'Estimated Zero Mean')

%estimatedExcerpt = zeros(length(excerpt), 1);
%estimatedExcerpt(mi_time-estimatedDelay) = excerpt(mi_time-estimatedDelay).*mi';

%clean = mixture;
%clean(delay:delay+length(excerpt) - 1) = clean(delay:delay+length(excerpt) - 1) - estimatedExcerpt;